import warnings
import sys
import copy
from pathlib import Path

import numpy as np
import pandas as pd

from skbio import diversity

from ecotools.util import guess_subsampling_level, elt_or_nothing

from ecotools.core.abundance_table import AbundanceTable
from ecotools.core.taxonomy_table import TaxonomyTable
from ecotools.core.metadata_table import MetadataTable
from ecotools.core.sequencing_data import SequencingData
from ecotools.core.phylogenetic_tree import PhylogeneticTree


class MetagenomicDS:
    def __init__(
            self,
            ds='project',
            outdir='.',
            abundance=None, metadata=None, taxonomy=None, fasta=None, tree=None,
            **kwargs
    ):
        self.outdir = Path(outdir)
        self.h5 = Path(outdir, ds+'.h5')
        self.figdir = Path(outdir, 'figures')
        
        self.abundance = AbundanceTable(data=abundance)
        self.taxonomy = TaxonomyTable(data=taxonomy)
        self.metadata = MetadataTable(data=metadata)
        self.sequences = SequencingData(data=fasta)
        self.tree = PhylogeneticTree(data=tree)

        self.distance_matrix = None

        self.coalesce()
        self.figdir.mkdir(parents=True, exist_ok=True)

    def __repr__(self):
        return '\n'.join(
            repr(data) for data in
            [self.abundance, self.taxonomy, self.metadata]
        )

    def __getattr__(self, key):

        if key in dir(self):
            return getattr(self, key)

        if not key.startswith('__') and 'abundance' in dir(self):
            abd_obj = getattr(self, 'abundance')

            try:
                return getattr(abd_obj, key)
            except AttributeError: 
                pass

        raise AttributeError(key)

    def copy(self):
        return copy.deepcopy(self)
    
    def n_samples(self):
        return self.abundance.data.shape[0]

    def n_otus(self):
        return self.abundance.data.shape[1]

    def factors(self):
        return self.metadata.qual_vars

    def covariates(self):
        return self.metadata.quant_vars

    def shapes(self):
        return {
            'otus': len(self.columns),
            'samples': len(self.index),
            'metadata': len(self.factors()) + len(self.covariates())
        }

    def get_column_format(self):
        abundance = self.abundance.get_column_format('group', 'OTU', 'value').reset_index()

        data = (abundance
                .merge(self.metadata.factor_data(), left_on='group', right_index=True)
                .merge(self.taxonomy.data, left_on='OTU', right_index=True)
                .set_index(['group', 'OTU']))

        return data
    
    def to_h5(self):
        if len(self.raw_sample_sizes) > len(self.index):
            warnings.warn('This method needs more work and might generate errors if an already processed dataset is used')
        if self.h5.is_file():
            self.h5.unlink()
        self.coalesce()
        self.abundance.to_h5(f'{self.outdir}/abundance.h5')
        self.taxonomy.to_h5(f'{self.outdir}/taxonomy.h5')
        self.metadata.to_h5(f'{self.outdir}/metadata.h5')

        if self.sequences is not None:
            self.sequences.to_h5(f'{self.outdir}/sequences.h5')

    def to_csv(self):
        self.coalesce()
        self.abundance.to_csv(f'{self.outdir}/abundance.csv')
        self.taxonomy.to_csv(f'{self.outdir}/taxonomy.csv')
        self.metadata.to_csv(f'{self.outdir}/metadata.csv')

    def coalesce(self, samples=True, otus=True, fasta=False):

        if samples:
            common_samples = np.intersect1d(self.abundance.index, self.metadata.index)
        
            if common_samples.size == 0:
                sys.exit("No common sample names between metadata and abundance table")

            self.abundance.subset_rows(common_samples)
            self.metadata.subset_rows(common_samples)

            if self.alpha_diversity is not None:
                self.alpha_diversity = self.alpha_diversity.loc[common_samples]
            if self.distance_matrix is not None:
                self.distance_matrix = self.distance_matrix.loc[common_samples, common_samples]
                

        if otus:
            is_common = self.abundance.columns.isin(self.taxonomy.index)
        
            if is_common.sum() == 0:
                sys.exit("Did not find any common OTUs between taxonomy and abundance table")

            common_otus = self.abundance.columns[is_common]
                
            self.abundance.subset_cols(common_otus)
            self.taxonomy.subset_rows(common_otus)
            # later: subset from the fasta + tree as well

        # Additional check: when abundance table was subsetted and we
        # - removed otus which created empty samples
        # - or removed samples which creates empty otus
        self.metadata.subset_rows(self.abundance.index)
        self.taxonomy.subset_rows(self.abundance.columns)

        if fasta:
            self.sequences.subset_otus(self.abundance.columns)
            self.tree.update_tree(self.sequences)
        
    def subset_samples(self, sample_names=None, sample_file=None):

        if sample_file is not None:
            sample_names = pd.read_csv(sample_file).iloc[:, 0]
        
        self.abundance.subset_rows(sample_names)
        self.coalesce(otus=False)

    def subset_otus(self, otus=None, taxa_files=None, taxa=None, clade=False):

        if taxa is not None:
            otus = self.taxonomy.get_ranks(taxa)

        if taxa_files is not None:
            if not isinstance(taxa_files, list):
                taxa_files = [taxa_files]

            ids = pd.concat([pd.read_csv(f).iloc[:, 0] for f in taxa_files])

            if ids.name.lower() == 'species':
                if clade:
                    hits = self.taxonomy.get_clade(ids)
                else:
                    hits = self.taxonomy.get_species(ids)
                otus = hits.index

            else:
                annotations = self.taxonomy.data[ids.name.title()]
                otus = self.columns[np.isin(annotations, ids)]

        self.abundance.subset_cols(otus)
        self.coalesce(samples=False)

    def group_samples(self, groups, fn='mean'):

        if isinstance(groups, str):
            groups = [groups]
            
        for i, group in enumerate(groups):
            if isinstance(group, str):
                groups[i] = self.metadata.data[group]

        group_names = [group.name for group in groups]

        self.abundance.data = self.abundance.data.groupby(groups).agg(fn)
        self.metadata.group_samples(group_names)

    def group_taxa(self, rank, discard_unknown=True):
        all_ranks = self.taxonomy.columns
        rank_idx = all_ranks.get_loc(rank) + 1
        ranks_to_keep = all_ranks[:rank_idx]
        rank_annot = self.taxonomy.data[rank]

        if discard_unknown:
            # Remove any rank with unknown label
            valid = ~(rank_annot.str.contains('uncultured|unclassified|unknown'))
            self.subset_otus(otus=self.columns[valid])
            
        self.abundance.data = self.abundance.data.groupby(rank_annot, axis=1, sort=False).agg(sum)
        tax = self.taxonomy.data.groupby(rank, sort=False)[ranks_to_keep].agg(elt_or_nothing)

        if tax.isnull().sum().sum() > 0:
            print('Some {}s have different upper level ranks'.format(rank))
            suspects = tax.index[tax.isnull().any(axis=1)]
            print(self.taxonomy.data.loc[self.taxonomy.data[rank].isin(suspects)])

            tax.dropna(inplace=True)
            self.abundance.data = self.abundance.data.loc[:, tax.index]

        self.taxonomy.data = tax

    def subsample(self, level=-1):
        if level < 0:
            sample_sums = self.abundance.data.sum(axis=1)
            level = guess_subsampling_level(sample_sums)
        self.abundance.subsample(level)        
        self.coalesce(otus=False)
        print('Subsampled at {}. {} samples remaining'.format(level, self.n_samples()))

    def compute_distance_matrix(self, metric='braycurtis', cache=False):

        sample_sums = self.abundance.data.sum(axis=1)
        if sample_sums.std() > 1:
            warnings.warn("Your samples are not of equal sizes. This is known to affect diversity calculation.", UserWarning)

        metric = metric.lower()
        print('Calculating {} distance for all samples ({})'.format(metric, self.data.shape[0]))

        output = Path(self.outdir, 'distances_{}_by-{}.npy'.format(metric, self.index.name))
        if output.is_file() and cache:
            dist = np.load(output)

        else:
            if metric in {'unweighted_unifrac', 'weighted_unifrac'}:
                self.sequences.subset_otus(self.abundance.columns)
                self.tree.update_tree(self.sequences)
                
                dist = diversity.beta_diversity(
                    metric, self.data, otu_ids=self.columns, tree=self.tree.data
                ).data
            else:
                dist = diversity.beta_diversity(metric, self.data).data
    
            np.save(output, dist)

        self.distance_matrix = pd.DataFrame(dist, index=self.index, columns=self.index)
        self.coalesce(otus=False)

    # def sort_otus(self):
    #     ordered_otus = self.abundance.data.sum().sort_values(ascending=False).index

    #     self.abundance.data = self.abundance.data[ordered_otus]
    #     self.taxonomy.data = self.taxonomy.data.loc[ordered_otus]
