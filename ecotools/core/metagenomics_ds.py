import copy
import inspect
import warnings
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from skbio import diversity

from ecotools.util import guess_subsampling_level, make_path
from ecotools.decorators import timer
from ecotools.rpy2_util import vegdist, to_phyloseq

from ecotools.core.abundance_table import AbundanceTable
from ecotools.core.taxonomy_table import TaxonomyTable
from ecotools.core.metadata_table import MetadataTable
from ecotools.core.sequencing_data import SequencingData
from ecotools.core.phylogenetic_tree import PhylogeneticTree


class MetagenomicDS:
    subclasses = {'abundance': AbundanceTable, 'taxonomy': TaxonomyTable,
                  'metadata': MetadataTable, 'sequences': SequencingData,
                  'tree': PhylogeneticTree}
    
    def __init__(
            self,
            ds='project',
            outdir='.', figdir=None, clade='',
            abundance=None, metadata=None, taxonomy=None, sequences=None, tree=None,
            distance_matrix=None, alpha_diversity=None, metric='bray', ordination='nmds'
    ):
        self.ds = ds
        self.outdir = Path(outdir)
        self.h5dir = Path(outdir, f'{ds}_h5-data')
        self.sample_axis = 'Group'
        self.feature_axis = 'OTU'

        if figdir is None:
            figdir = Path(outdir, 'figures')
        self.figdir = Path(figdir)
            
        self.run = {
            'metric': metric,
            'ordination': ordination,
            'clade': clade,
            'factors': None,
            'strata': None
        }

        self.io = {
            'outdir': Path(outdir), 'h5dir': Path(outdir, f'{ds}_h5-data'), 'figdir': figdir,
            'distances': lambda args: make_path(outdir, '.csv', 'distances', metric, *args),
            'ordination': lambda args: make_path(outdir, '.csv', ordination, metric, *args)
        }
        
        # Set special attributes (sub-classes).
        # If the user provides the data, calls the constructor,
        # otherwise set the attribute directly
        frame = inspect.currentframe()
        args, _, _, locals_ = inspect.getargvalues(frame)
        
        for arg in MetagenomicDS.subclasses:
            class_obj = MetagenomicDS.subclasses[arg]
            if isinstance(locals_[arg], class_obj):
                setattr(self, arg, locals_[arg])
            else:
                setattr(self, arg, class_obj(data=locals_[arg]))

        self.distance_matrix = distance_matrix

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
        kwargs = {attr: copy.copy(getattr(self, attr))
                  for attr in inspect.signature(self.__init__).parameters.keys()
                  if hasattr(self, attr)}

        obj = MetagenomicDS(**kwargs)
        obj.io = self.io
        obj.run = self.run

        return obj

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

    def get_column_format(self, tax=True, meta=True):
        data = self.abundance.get_column_format('group', 'OTU', 'value').reset_index()

        if meta:
            data = data.merge(self.metadata.factor_data(), left_on='group', right_index=True)
        if tax:
            data = data.merge(self.taxonomy.data, left_on='OTU', right_index=True)

        data.set_index(['group', 'OTU'], inplace=True)

        return data
    
    def to_h5(self):
        self.h5dir.mkdir(exist_ok=True)
        
        if len(self.raw_sample_sizes) > len(self.index):
            warnings.warn('This method needs more work and might generate errors if an already processed dataset is used')
        self.coalesce()
        
        self.abundance.to_h5(outdir=self.io["h5dir"], filename='abundance.h5')
        self.taxonomy.to_h5(outdir=self.io["h5dir"], filename='taxonomy.h5')
        self.metadata.to_h5(outdir=self.io["h5dir"], filename='metadata.h5')

        if self.sequences is not None:
            self.sequences.to_h5(outdir=self.h5dir, filename='sequences.h5')
        if self.tree is not None:
            self.tree.to_h5(outdir=self.h5dir, filename='tree.h5')

    def to_csv(self, outdir=None):
        if outdir is None:
            outdir = self.io['outdir']
        self.coalesce()
        self.abundance.to_csv(f'{outdir}/abundance.csv')
        self.taxonomy.to_csv(f'{outdir}/taxonomy.csv')
        self.metadata.to_csv(f'{outdir}/metadata.csv')

    def to_phyloseq(self):
        try:
            tree_path = self.tree.tree_path
        except AttributeError:
            tree_path = None
        phylo_obj = to_phyloseq(self.abundance.data,
                                self.taxonomy.data,
                                self.metadata.data,
                                tree_path=tree_path)
        return phylo_obj
        
    def coalesce(self, samples=True, otus=True, fasta=False):

        if samples:
            common_samples = np.intersect1d(self.abundance.index,
                                            self.metadata.index)
        
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
        
    def subset_samples(self, sample_names=None, sample_file=None, inplace=True, **kwargs):

        if sample_file is not None:
            sample_names = pd.read_csv(sample_file).iloc[:, 0]

        if kwargs:
            meta_info = self.metadata.factor_data(kwargs.keys())
            sample_cond = (meta_info == list(kwargs.values())).all(axis=1)
            sample_names = self.abundance.index[sample_cond]
            
        if not inplace:
            mg = self.copy()
            mg.abundance.subset_rows(sample_names)
            mg.coalesce()
            return mg
            
        self.abundance.subset_rows(sample_names)
        self.coalesce(otus=False)

    def subset_otus(self, otus=None, taxa=None, clade=False, sep='\t', header='infer'):

        if taxa is not None:
            otus = self.taxonomy.get_otus_from_pairs(taxa)

        if isinstance(otus, Path) or (isinstance(otus, list) and isinstance(otus[0], Path)):
            if not isinstance(otus, list):
                otus = [otus]

            ids = pd.concat([pd.read_csv(f, sep=sep, header=header).iloc[:, 0] for f in otus])

            if ids.name.title() == 'Species':
                if clade:
                    hits = self.taxonomy.get_clade(ids)
                else:
                    hits = self.taxonomy.get_species(ids)
                otus = hits.index

            elif ids.name.title() in TaxonomyTable.ranks:
                annotations = self.taxonomy.data[ids.name.title()]
                otus = self.columns[np.isin(annotations, ids)]
            else:
                otus = ids.values

        self.abundance.subset_cols(otus)
        self.coalesce(samples=False)

    def group_samples(self, *groups, fn='mean'):

        groups = list(groups)
        for i, group in enumerate(groups):
            if isinstance(group, str):
                groups[i] = self.metadata.data[group]

        self.abundance.group_samples(groups, fn=fn)
        self.metadata.group_samples([group.name for group in groups], fn=fn)
        self.sample_axis = '-'.join(g.name for g in groups)

    def group_taxa(self, rank, **kwargs):
        mapping = self.taxonomy.data[rank.title()]
        self.taxonomy.group_taxa(rank, **kwargs)
        self.abundance.data = (
            self.abundance.data.groupby(mapping, axis=1).agg(sum)
        )
        self.coalesce(samples=False)
        self.feature_axis = rank

    def subsample(self, level=-1, plot=False):
        if level < 0:
            sample_sums = self.abundance.data.sum(axis=1)
            level = guess_subsampling_level(sample_sums, plot=plot)
        n_before = self.n_samples()

        self.abundance.subsample(level)        
        self.coalesce(otus=False)
        print('Subsampled at {}. {}/{} samples remaining'
              .format(level, self.n_samples(), n_before))
        return level

    @timer
    def compute_distance_matrix(self, cache=False, vegan=True, strata=None):

        sample_sums = self.abundance.data.sum(axis=1)
        if sample_sums.std() > 1:
            warnings.warn("Your samples are not of equal sizes. This is known to affect diversity calculation.", UserWarning)

        output = self.io['distances']((self.run['clade'], strata))
    
        if output.is_file() and cache:
            self.distance_matrix = pd.read_csv(output, index_col=[0, 1])[self.run['metric']]

        if vegan and not self.run['metric'].endswith('unifrac'):
            # Use R vegan package
            self.distance_matrix = vegdist(self.data, self.run['metric'], r_obj=False)
        else:
            # Use skbio
            if self.run['metric'].lower().endswith('unifrac'):
                if self.tree.data is None:
                    if self.tree_path is None:
                        self.tree.compute(sequences=self.sequences)
                    else:
                        self.tree.load()
                self.tree.root()

                dist = diversity.beta_diversity(
                    self.run['metric'], self.abundance.data,
                    otu_ids=self.columns, tree=self.tree.data
                )
            else:
                dist = diversity.beta_diversity(self.run['metric'], self.data)

            self.distance_matrix = pd.DataFrame(dist.data, index=self.index, columns=self.index).stack()

        self.distance_matrix.to_csv(output)

