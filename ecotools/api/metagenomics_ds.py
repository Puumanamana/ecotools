import sys
import copy
from pathlib import Path

import numpy as np
import pandas as pd

from ecotools.api.base_data_classes import AbundanceTable, TaxonomyTable, MetadataTable

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
            fasta_path='',
            meta_path='', qual_vars=None
    ):
        self.outdir = Path(outdir)
        self.h5 = Path(outdir, ds+'.h5')
        self.figdir = Path(outdir, 'figures')
        
        self.abundance = AbundanceTable(path=abd_path, ds=ds, outdir=outdir)
        self.taxonomy = TaxonomyTable(path=tax_path, species_path=species_path, ds=ds, outdir=outdir)
        self.metadata = MetadataTable(*qual_vars, path=meta_path, ds=ds, outdir=outdir)
        self.unify()
        self.figdir.mkdir(parents=True, exist_ok=True)

    def __repr__(self):
        return '\n'.join(
            data.__repr__() for data in
            [self.abundance, self.taxonomy, self.metadata]
        )

    def copy(self):
        return copy.deepcopy(self)

    def otus(self):
        return self.abundance.data.columns

    def samples(self):
        return self.abundance.data.index

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
    
    def to_h5(self):
        self.unify()
        self.abundance.to_h5()
        self.taxonomy.to_h5()
        self.metadata.to_h5()

    def to_csv(self):
        self.unify()
        self.abundance.to_csv(Path(self.outdir, "abundance.csv"))
        self.taxonomy.to_csv(Path(self.outdir, "taxonomy.csv"))
        self.metadata.to_csv(Path(self.outdir, "metadata.csv"))

    def unify(self):
        common_samples = np.intersect1d(self.abundance.data.index, self.metadata.data.index)

        if common_samples.size == 0:
            sys.exit("Did not find any common sample names between metadata and abundance table")

        self.subset_samples(sample_names=common_samples)

        common_otus = np.intersect1d(self.abundance.data.columns, self.taxonomy.data.index)
        
        if common_otus.size == 0:
            sys.exit("Did not find any common OTUs between taxonomy and abundance table")

        self.subset_otus(otus=common_otus)

    def subset_samples(self, sample_names=None, sample_file=None):

        if sample_file is not None:
            sample_names = pd.read_csv(sample_file).iloc[:, 0]
        
        self.abundance.data = self.abundance.data.loc[sample_names]
        self.metadata.data = self.metadata.data.loc[sample_names]

        null_otus = self.abundance.data.sum() == 0

        if sum(null_otus) > 0:
            self.abundance.data = self.abundance.data.loc[:, ~null_otus]
            self.taxonomy.data = self.taxonomy.data.loc[~null_otus]

    def subset_otus(self, otus=None, taxa_file=None):

        if taxa_file is not None:
            ids = pd.read_csv(taxa_file).iloc[:, 0]

            if ids.name.lower() == 'species':
                hits = self.taxonomy.subset_species(ids)
                otus = hits.index

            else:
                annotations = self.taxonomy.data[ids.name.title()]
                otus = self.otus()[np.isin(annotations, ids)]

        self.abundance.data = self.abundance.data.loc[:, otus]
        self.taxonomy.data = self.taxonomy.data.loc[otus]

        null_samples = self.abundance.data.sum(axis=1) == 0

        if sum(null_samples) > 0:
            self.abundance.data = self.abundance.data.loc[~null_samples]
            self.metadata.data = self.metadata.data.loc[~null_samples]

    def group_samples(self, columns, fn):
        groups = self.metadata.data[columns]

        self.abundance.data = self.abundance.data.groupby(groups).agg(fn)

        self.metadata.group_samples(columns)

    def group_taxa(self, rank):
        all_ranks = self.taxonomy.data.columns
        rank_idx = all_ranks.get_loc(rank) + 1
        ranks_to_keep = all_ranks[:rank_idx]

        rank_annot = self.taxonomy.data[rank]

        valid = ~(rank_annot.str.contains('uncultured|unclassified|unknown'))

        self.subset_otus(otus=self.otus()[valid])
        self.abundance.data = self.abundance.data.groupby(rank_annot, axis=1).agg(sum)
        self.taxonomy.data = self.taxonomy.data.groupby(rank)[ranks_to_keep].agg(elt_or_nothing)

        if self.taxonomy.data[rank].isnull().sum() > 0:
            print('Some {}s have different upper level ranks'.format(rank))
            import ipdb;ipdb.set_trace()

    def get_top_otus(self, n=-1):
        sorted_otus = self.abundance.data.sum().sort_values(ascending=False).index
        
        if n <=0:
            return sorted_otus

        return sorted_otus[:n]

    
    def preprocess(self, factor=None, taxa_file=None, norm=False, rank=None, top=-1):

        if taxa_file is not None:
            self.subset_otus(taxa_file=taxa_file)

        if rank is not None:
            # Only show annotation at the `rank` level
            self.group_taxa(rank)

        if norm:
            # Normalize by sample sum
            self.abundance.normalize()

        if factor is not None:
            # Mean of OTU abundance in each group
            self.group_samples(factor, 'mean')

        # Only show the most `top` abundant taxa. -1 is just sorting the taxa
        top_otus = self.get_top_otus(top)
        top = len(top_otus) # in case top=-1
        self.subset_otus(otus=top_otus)
        
