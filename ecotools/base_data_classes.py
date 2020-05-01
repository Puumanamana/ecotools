import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import h5py

from skbio import diversity
from skbio.stats import subsample_counts

from Bio.SeqIO.FastaIO import SimpleFastaParser

def get_header(path, sep='\t'):
    return next(open(path)).split(sep)

def get_str_dtype(data):

    def str_len(x): return len(str(x))
    
    if len(data.shape) == 1:
        max_len = data.map(str_len).max()
    else:
        max_len = data.applymap(str_len).max().max()

    return 'S{}'.format(max_len)

class Data:
    def __init__(self, ds='project', path='', outdir='.'):
        self.ds = ds
        self.path = Path(path)
        self.outdir = outdir
        self.h5 = Path(outdir, ds+'.h5')
        self.data = None

    def valid_input(self):
        if not self.path.is_file():
            msg = 'Could not find {}'.format(self.path.name)
            warnings.warn(msg, UserWarning)
            return False
        return True

    def subset_rows(self, x):
        self.data = self.data.loc[x]

    def subset_cols(self, x):
        self.data = self.data.loc[:, x]        

    def to_csv(self, path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.data.to_csv(path)

class MetadataTable(Data):
    def __init__(self, *qual_vars, **kwargs):
        Data.__init__(self, **kwargs)
        self.quant_vars = None
        self.qual_vars = None
        
        self.load(dtypes={v: str for v in qual_vars})

    def __repr__(self):
        return "Metadata table: {} samples, {} columns ({} factors, {} covariates)".format(
            *self.data.shape, len(self.qual_vars), len(self.quant_vars)
        )

    def load(self, dtypes):
        if self.h5.is_file():
            self.load_h5()
            return

        if not self.valid_input():
            return

        if self.path.suffix not in {'.csv', '.tsv'}:
            warnings.warn('Unknown extension for file: {}. Choices are ".csv" or ".tsv"'
                          .format(self.path.name), UserWarning)

        sep = {'.csv': ',', '.tsv': '\t'}.get(self.path.suffix, '\t')

        self.data = pd.read_csv(self.path, index_col=0, sep=sep, dtype=dtypes)
        self.data.index.name = 'SampleID'

        self.quant_vars = self.data.select_dtypes(include='number').columns.to_numpy()
        self.qual_vars = self.data.drop(self.quant_vars, axis=1).columns.to_numpy()

        if not self.valid_numeric_dtypes():
            sys.exit('There was a problem when parsing metadata. Some numeric columns contain invalid "numbers"')

    def valid_numeric_dtypes(self):
        numeric_columns = self.data.select_dtypes(include='number').columns
        not_numeric_quant_vars = [col for col in self.quant_vars if col not in numeric_columns]

        if not_numeric_quant_vars:
            return False
        return True
            
    def load_h5(self): 
        with h5py.File(self.h5, 'r') as handle:
            table_qual = handle.get('metadata_qual')[:].astype(str)
            table_quant = handle.get('metadata_quant')[:]
            samples = handle.get('samples')[:].astype(str)
            self.quant_vars = handle.get('quant_vars')[:].astype(str)
            self.qual_vars = handle.get('qual_vars')[:].astype(str)

        self.data = pd.concat([
            pd.DataFrame(table_qual, index=samples, columns=self.qual_vars),
            pd.DataFrame(table_quant, index=samples, columns=self.quant_vars),
        ], axis=1)
    
    def to_h5(self):
        qual_dtype = get_str_dtype(self.data[self.qual_vars])
        col_dtype = get_str_dtype(self.data.columns)
        sample_dtype = get_str_dtype(self.data.index)

        with h5py.File(self.h5, 'a') as handle:
            handle.create_dataset(
                'metadata_qual',
                data=self.data[self.qual_vars].to_numpy().astype(qual_dtype)
            )
            handle.create_dataset(
                'metadata_quant',
                data=self.data[self.quant_vars].to_numpy()
            )            
            handle.create_dataset(
                'qual_vars',
                data=self.qual_vars.astype(col_dtype)
            )
            handle.create_dataset(
                'quant_vars',
                data=self.quant_vars.astype(col_dtype)
            )

            if not 'samples' in handle:
                handle.create_dataset(
                    'samples',
                    data=self.data.index.to_numpy().astype(sample_dtype)
                )

    def group_samples(self, factors):
        # keep the columns with a unique value when aggregated
        cols_to_drop = self.qual_vars[self.data
                        .groupby(factors)[self.qual_vars]
                        .agg(lambda x: len(set(x)))
                        .max() > 1]

        print('Dropping {}: multiple values when grouping'.format(
            ','.join(cols_to_drop.tolist()))
        )
        self.quant_vars = np.setdiff1d(self.quant_vars, cols_to_drop)
        self.qual_vars = np.setdiff1d(self.qual_vars, cols_to_drop)

        self.data = self.data.groupby(factors).agg('first').drop(cols_to_drop, axis=1)

    def factor_data(self, col=None):
        if col is None:
            return self.data[self.qual_vars]
        return self.data[col]

    def subset_rows(self, x):
        self.data = self.data.loc[x]
        self.data = self.data.loc[:, self.data.count() > 0]

    def add_var(self, name, values):
        try:
            values = pd.to_numeric(values)
            self.quant_vars = np.append(self.quant_vars, name)
        except:
            self.qual_vars = np.append(self.qual_vars, name)

        self.data.loc[:, name] = values

    def drop_vars(self, *names):
        names = list(names)
        self.quant_vars = np.setdiff1d(self.quant_vars, names)
        self.qual_vars = np.setdiff1d(self.qual_vars, names)
        self.data.drop(names, axis=1, inplace=True)

class AbundanceTable(Data):

    def __init__(self, **kwargs):
        Data.__init__(self, **kwargs)
        self.load()
        self.raw_sample_sizes = self.data.sum(axis=1)
        self.alpha_diversity = None
        self.beta_diversity = None

    def __repr__(self):
        return "Abundance table: {} samples, {} OTUs".format(*self.data.shape)

    def load(self):
        if self.h5.is_file():
            self.load_h5()
            return

        if not self.valid_input():
            return
        
        header = get_header(self.path)
        dtypes = dict((col, str) if col in ['Group'] else (col, int) for col in header)

        self.data = pd.read_csv(self.path, index_col='Group', sep='\t', dtype=dtypes,
                                keep_default_na=False, low_memory=False,
                                usecols=lambda x: x not in ['label', 'numOtus'])
        self.data.index.name = 'SampleID'

    def to_h5(self):
        sample_dtype = get_str_dtype(self.data.index)
        otu_dtype = get_str_dtype(self.data.columns)
        
        with h5py.File(self.h5, 'a') as handle:
            handle.create_dataset(
                'abundance',
                data=self.data.to_numpy()
            )

            if not 'samples' in handle:
                handle.create_dataset(
                    'samples',
                    data=self.data.index.to_numpy().astype(sample_dtype)
                )

            if not 'OTUs' in handle:
                handle.create_dataset(
                    'OTUs',
                    data=self.data.columns.to_numpy().astype(otu_dtype)
                )

    def load_h5(self):
        with h5py.File(self.h5, 'r') as handle:
            table = handle.get('abundance')[:]
            samples = handle.get('samples')[:].astype(str)
            OTUs = handle.get('OTUs')[:].astype(str)

        self.data = pd.DataFrame(table, index=samples, columns=OTUs)

    def sample_sizes(self):
        return self.data.sum(axis=1)

    def subset_rows(self, x):
        self.data = self.data.loc[x]
        self.data = self.data.loc[:, self.data.sum() > 0]
    
    def subset_cols(self, x):
        self.data = self.data.loc[:, x]
        self.data = self.data.loc[self.data.sum(axis=1) > 0, :]
    
    def normalize(self):
        self.data = ((self.data.T) / self.raw_sample_sizes.loc[self.data.index]).T

    def wisconsin(self):
        self.data = self.data / self.data.max(axis=0)
        self.data = (self.data.T / self.data.sum(axis=1)).T

    def subsample(self, level):
        dropped = []

        for (i, row) in enumerate(self.data.to_numpy()):
            try:
                row_subsampled = subsample_counts(row, level)
            except ValueError:
                dropped.append(i)
                continue
                
            self.data.iloc[i] = row_subsampled

        self.data.drop(self.data.index[dropped], inplace=True)

    def get_most_abundant_otus(self, thresh=0.01):
        proportions = (self.data.T / self.data.sum(axis=1)).max(axis=1)
        main_otus = proportions.index[proportions > thresh]

        return main_otus

    def compute_alpha_diversity(self):
        alpha_div = pd.DataFrame({
            'n_otus': self.data.sum(axis=1),
            'richness': (self.data > 0).sum(axis=1),
            'chao': self.data.apply(diversity.alpha.chao1, axis=1),
            'shannon': self.data.apply(diversity.alpha.shannon, axis=1),
        })

        self.alpha_diversity = alpha_div

    def compute_distance_matrix(self, tree=None, metric='braycurtis', cache=False, force_subsampling=True):

        sample_sums = self.data.sum(axis=1)
        if sample_sums.std() > 1:
            warnings.warn("Your samples are not of equal sizes. This is known to affect diversity calculation. Performing subsampling. You can choose to not subsample at your own risks by setting force_subsampling to False.", UserWarning)

            if force_subsampling:
                subsampling_level = max(sample_sums.quantile(0.1, axis=1), 5000)
                self.subsample(subsampling_level)
        
        metric = metric.lower()
        print('Calculating {} distance for all samples ({})'.format(metric, self.data.shape[0]))

        output = Path(self.outdir, 'distances_{}_by-{}.npy'.format(metric, self.data.index.name))
        if output.is_file() and cache:
            dist = np.load(output)

        else:
            if metric in {'unweighted_unifrac', 'weighted_unifrac'}:
                dist = diversity.beta_diversity(
                    metric, self.data, otu_ids=self.data.columns, tree=tree
                ).data
            else:
                dist = diversity.beta_diversity(metric, self.data).data
                
            np.save(output, dist)

        self.beta_diversity = pd.DataFrame(dist, index=self.data.index, columns=self.data.index)

class TaxonomyTable(Data):
    ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    
    def __init__(self, species_path='', **kwargs):
        Data.__init__(self, **kwargs)
        self.species_path = species_path
        self.load()

    def __repr__(self):
        return "Taxonomy table: {} OTUs, {} levels".format(*self.data.shape)

    def load(self):
        if self.h5.is_file():
            self.load_h5()
            return

        if not self.valid_input():
            return

        self.data = (pd.read_csv(self.path, index_col=0, sep='\t').Taxonomy
                     .str.strip(';')
                     .str.split(';', expand=True))

        if self.species_path:
            species_info = pd.read_csv(self.species_path, index_col=0).Species
            self.data = pd.concat([self.data, species_info], axis=1)

        self.data.columns = TaxonomyTable.ranks[:self.data.shape[1]]

    def load_h5(self, outdir='.', ds='project'):
        with h5py.File(self.h5, 'r') as handle:
            table = handle.get('taxonomy')[:]
            levels = handle.get('tax_levels')[:].astype(str)
            OTUs = handle.get('OTUs')[:].astype(str)

        self.data = pd.DataFrame(table, index=OTUs, columns=levels)
    
    def to_h5(self):
        dtype = get_str_dtype(self.data)
        otu_dtype = get_str_dtype(self.data.index)

        with h5py.File(self.h5, 'a') as handle:
            handle.create_dataset(
                'taxonomy',
                data=self.data.to_numpy().astype(dtype)
            )
            handle.create_dataset(
                'tax_levels',
                data=self.data.columns.to_numpy().astype('S8')
            )

            if not 'OTUs' in handle:
                handle.create_dataset(
                    'OTUs',
                    data=self.data.index.to_numpy().astype(otu_dtype)
                )

    def get_ranks(self, info):
        '''
        Get either of the rank, values in info (logical or)
        '''

        info = [(rank.title(), pd.Series(vals).str.lower()) for (rank, vals) in info]
        data_lower = self.data.fillna('').applymap(lambda x: x.lower())
        
        conditions = [
            data_lower[rank].isin(vals)
            for (rank, vals) in info
        ]

        if len(conditions) == 1:
            condition = conditions[0]
        else:
            condition = np.logical_or(*conditions)

        return self.data.index[condition]

    def get_clade(self, other, threshold=0.5):
        if not self.species_path:
            sys.exit('Error: No species information')

        species = self.data.Species.dropna().str.split('/').explode()
        species = self.data.reindex(index=species.index).Genus + ' ' + species

        is_clade = (species.str.lower()
                    .isin(other.str.lower())
                    .groupby(level=0).agg('mean'))
        
        clade_otus = is_clade.index[is_clade > threshold]
        
        hits = self.data.Species.loc[clade_otus]

        return hits

    def get_species(self, other):
        return self.get_clade(other, threshold=0)

    
class DNAsequences(Data):

    def __init__(self, **kwargs):
        Data.__init__(self, **kwargs)
        self.load()

    def __repr__(self):
        return "OTUs: {} sequences".format(len(self.data))

    def load_h5(self, outdir='.', ds='project'):
        print('Not implemented')
    
    def load(self):
        if self.h5.is_file():
            self.load_h5()
            return

        if not self.valid_input():
            return
        
        with h5py.File(self.path, 'r') as in_handle:
            self.data = []
            for title, seq in SimpleFastaParser(in_handle):
                (name, descr) = title.split()
                self.data.append({'name': name, 'descr': descr, 'seq': seq})

# class PhyloTree(Data):

#     def __init__(self, **kwargs):
#         Data.__init__(self, **kwargs)
#         self.load()

#     def __repr__(self):
#         return repr(self.data)

#     def load_h5(self, outdir='.', ds='project'):
#         print('Not implemented')
    
#     def load(self):
#         if self.h5.is_file():
#             self.load_h5()
#             return

#         if not self.valid_input():
#             return
        
#         self.data = TreeNode.read(str(self.path))
