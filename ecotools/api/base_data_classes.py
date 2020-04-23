import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
from skbio import diversity

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
            warnings.warn('Could not find {}'.format(self.path.name), UserWarning)
            return False
        return True

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

    def add_var(self, name, values):
        try:
            values = pd.to_numeric(values)
            self.quant_vars = np.append(self.quant_vars, name)
        except:
            self.qual_vars = np.append(self.qual_vars, name)

        self.data.loc[:, name] = values

class AbundanceTable(Data):

    def __init__(self, **kwargs):
        Data.__init__(self, **kwargs)
        self.load()
        self.raw_sample_sizes = self.data.sum(axis=1)
        self.alpha_diversity = None

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

    def normalize(self):
        self.data = ((self.data.T) / self.raw_sample_sizes.loc[self.data.index]).T

    def get_most_abundant_otus(self, thresh=0.01):
        proportions = (self.data / self.data.sum()).min()
        main_otus = proportions.index[proportions > thresh]

        return main_otus

    def calc_diversity(self):

        alpha_div = pd.DataFrame({
            'n_otus': self.data.sum(axis=1),
            'richness': (self.data > 0).sum(axis=1),
            'chao': self.data.apply(diversity.alpha.chao1, axis=1),
            'shannon': self.data.apply(diversity.alpha.shannon, axis=1),
        })

        self.alpha_diversity = alpha_div

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

    def subset_species(self, other, threshold=0.5):
        if not self.species_path:
            sys.exit('Error: No species information')

        species = self.data.Species.dropna().str.split('/').explode()
        species = self.data.reindex(index=species.index).Genus + ' ' + species

        is_pathogenic = (species.str.lower()
                         .isin(other.str.lower())
                         .groupby(level=0).agg('mean'))
        
        patho_otus = is_pathogenic.index[is_pathogenic >= threshold]
        
        hits = self.data.Species.loc[patho_otus]

        return hits

class DNAsequences(Data):

    def __init__(self, path):
        self.path = Path(path)
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
