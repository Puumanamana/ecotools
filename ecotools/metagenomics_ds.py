import warnings
import sys
import copy
from pathlib import Path

import numpy as np
import pandas as pd

from skbio import diversity

from ecotools.base_data_classes import AbundanceTable, TaxonomyTable, MetadataTable, SequenceData
from ecotools.plotting.heatmap import clustermap
from ecotools.util import guess_subsampling_level

def elt_or_nothing(l):
    if len(set(l)) == 1:
        return l[0]
    return None

class MetagenomicDS:
    def __init__(
            self,
            ds='project',
            outdir='.',
            abd_path='',
            tax_path='', species_path='',
            fasta_path=None,
            tree_path=None,
            meta_path='', qual_vars=None
    ):
        self.outdir = Path(outdir)
        self.h5 = Path(outdir, ds+'.h5')
        self.figdir = Path(outdir, 'figures')
        
        self.abundance = AbundanceTable(path=abd_path, ds=ds, outdir=outdir)
        self.taxonomy = TaxonomyTable(path=tax_path, species_path=species_path, ds=ds, outdir=outdir)
        self.seq_data = SequenceData(fasta_path=fasta_path, tree_path=tree_path, outdir=outdir)
        self.metadata = MetadataTable(*qual_vars, path=meta_path, ds=ds, outdir=outdir)
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

        if 'abundance' in dir(self):
            abd_obj = getattr(self, 'abundance')

            if key in dir(abd_obj):
                return getattr(abd_obj, key)

        if not key.startswith('__'):
            print('{} not found in {}'.format(key, self.__class__))

        raise AttributeError

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
            'otus': len(self.otus()),
            'samples': len(self.samples()),
            'metadata': len(self.factors()) + len(self.covariates())
        }

    def get_column_format(self):
        abundance = self.abundance.get_column_format('group', 'OTU', 'abundance').reset_index()

        data = (abundance
                .merge(self.metadata.data, left_on='group', right_index=True)
                .merge(self.taxonomy.data, left_on='OTU', right_index=True)
                .set_index(['group', 'OTU']))

        return data
    
    def to_h5(self):
        if len(self.raw_sample_sizes) > len(self.samples()):
            warnings.warn('This method needs more work and might generate errors if an already processed dataset is used')
        if self.h5.is_file():
            self.h5.unlink()
        self.coalesce()
        self.abundance.to_h5()
        self.taxonomy.to_h5()
        self.metadata.to_h5()

    def to_csv(self):
        self.coalesce()
        self.abundance.to_csv(Path(self.outdir, "abundance.csv"))
        self.taxonomy.to_csv(Path(self.outdir, "taxonomy.csv"))
        self.metadata.to_csv(Path(self.outdir, "metadata.csv"))

    def coalesce(self, samples=True, otus=True):

        if samples:
            common_samples = np.intersect1d(self.abundance.data.index, self.metadata.data.index)
        
            if common_samples.size == 0:
                sys.exit("No common sample names between metadata and abundance table")

            self.abundance.subset_rows(common_samples)
            self.metadata.subset_rows(common_samples)

            if self.alpha_diversity is not None:
                self.alpha_diversity = self.alpha_diversity.loc[common_samples]
            if self.distance_matrix is not None:
                self.distance_matrix = self.distance_matrix.loc[common_samples, common_samples]
                

        if otus:
            common_otus = np.intersect1d(self.abundance.data.columns, self.taxonomy.data.index)
        
            if common_otus.size == 0:
                sys.exit("Did not find any common OTUs between taxonomy and abundance table")

            self.abundance.subset_cols(common_otus)
            self.taxonomy.subset_rows(common_otus)
            # later: subset from the fasta + tree as well

        # Additional check: when abundance table was subsetted and we
        # - removed otus which created empty samples
        # - or removed samples which creates empty otus
        self.metadata.subset_rows(self.abundance.data.index)
        self.taxonomy.subset_rows(self.abundance.data.columns)
        
    def subset_samples(self, sample_names=None, sample_file=None):

        if sample_file is not None:
            sample_names = pd.read_csv(sample_file).iloc[:, 0]
        
        self.abundance.subset_rows(sample_names)
        self.coalesce(otus=False)

    def subset_otus(self, otus=None, taxa_files=None, taxa=None, clade=False):

        if taxa is not None:
            otus = self.taxonomy.get_ranks(taxa)

        if taxa_files is not None:
            if type(taxa_files) == str:
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
                otus = self.otus()[np.isin(annotations, ids)]

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
        all_ranks = self.taxonomy.data.columns
        rank_idx = all_ranks.get_loc(rank) + 1
        ranks_to_keep = all_ranks[:rank_idx]
        rank_annot = self.taxonomy.data[rank]

        if discard_unknown:
            # Remove any rank with unknown label
            valid = ~(rank_annot.str.contains('uncultured|unclassified|unknown'))
            self.subset_otus(otus=self.otus()[valid])
            
        self.abundance.data = self.abundance.data.groupby(rank_annot, axis=1).agg(sum)
        self.taxonomy.data = self.taxonomy.data.groupby(rank)[ranks_to_keep].agg(elt_or_nothing)

        if self.taxonomy.data[rank].isnull().sum() > 0:
            print('Some {}s have different upper level ranks'.format(rank))
            import ipdb;ipdb.set_trace()

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

        output = Path(self.outdir, 'distances_{}_by-{}.npy'.format(metric, self.data.index.name))
        if output.is_file() and cache:
            dist = np.load(output)

        else:
            if metric in {'unweighted_unifrac', 'weighted_unifrac'}:
                self.seq_data.update_tree(self.otus())
                
                dist = diversity.beta_diversity(
                    metric, self.data, otu_ids=self.data.columns, tree=self.seq_data.load_tree()
                ).data
            else:
                dist = diversity.beta_diversity(metric, self.data).data
    
            np.save(output, dist)

        self.distance_matrix = pd.DataFrame(dist, index=self.samples(), columns=self.samples())
        self.coalesce(otus=False)

    def preprocess(self, relabund=False, rank=None, top=-1,
                   taxa_files=None, taxa=None, clade=False):

        if taxa_files is not None or taxa is not None:
            self.subset_otus(taxa_files=taxa_files, clade=clade, taxa=taxa)

        if rank is not None:
            # Only show annotation at the `rank` level
            self.group_taxa(rank)

        if relabund:
            # Normalize by sample sum
            self.abundance.to_relative_abundance()

    def sort_otus(self):
        ordered_otus = self.abundance.data.sum().sort_values(ascending=False).index

        self.abundance.data = self.abundance.data[ordered_otus]
        self.taxonomy.data = self.taxonomy.data.loc[ordered_otus]

    def clustermap(self, rank=None, top=None, **kwargs):
        mg = self.copy()
        mg.to_relative_abundance()
        mg.taxonomy.clean_labels(trim=True)

        if rank is not None:
            mg.group_taxa(rank)

        if top is None and mg.n_otus() > 100:
            top = 100

        if top is not None:
            print('Getting top 100 most prevalent otus')
            otus = mg.select_otus(criteria='prevalence', n=100)
            mg.subset_otus(otus=otus)

        if 'output' in kwargs:
            kwargs['output'] = str(kwargs['output']).replace('.html', '_top-{}.html'.format(top))

        clustermap(mg, standardize=True, cluster_rows=True, cluster_cols=True, **kwargs)

        
