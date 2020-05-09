import sys
import warnings

from pathlib import Path
import subprocess
from tempfile import mkdtemp

import numpy as np
import pandas as pd
import h5py

import sklearn.preprocessing
from skbio import diversity
from skbio.stats import subsample_counts
from skbio.tree import TreeNode

from Bio import SeqIO, AlignIO
from Bio.Phylo.Applications import PhymlCommandline

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

    def get_column_format(self, index_name, col_name, value_name):
        data = self.data.stack().rename(value_name)
        data.index.names = [index_name, col_name]
        return data

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
        qual_vars = np.setdiff1d(self.qual_vars, factors)
        
        cols_to_drop = qual_vars[self.data
                                 .groupby(factors)[qual_vars]
                                 .agg(lambda x: len(set(x)))
                                 .max() > 1]

        self.qual_vars = np.setdiff1d(self.qual_vars, cols_to_drop)

        agg = dict((col, 'mean') if col in self.quant_vars
                   else (col, 'first') for col in self.data.columns)
        
        self.data = self.data.groupby(factors).agg(agg).drop(cols_to_drop, axis=1)
        self.data.index.names = factors

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
        self.distance_matrix = None

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

    def otus(self):
        return self.data.columns

    def samples(self):
        return self.data.index
    
    def sample_sizes(self):
        return self.data.sum(axis=1)

    def subset_rows(self, x):
        self.data = self.data.loc[x]
        self.data = self.data.loc[:, self.data.sum() > 0]
    
    def subset_cols(self, x):
        self.data = self.data.loc[:, x]
        self.data = self.data.loc[self.data.sum(axis=1) > 0, :]

    def to_relative_abundance(self):
        if self.data.shape[0] < len(self.raw_sample_sizes):
            print('Cannot normalize by sample size: data already aggregated')
            return
        self.data = (self.data.T / self.raw_sample_sizes).T
        
    def normalize(self, method, axis='features', inplace=True):
        data = self.data
        if axis == 'samples':
            data = data.T

        if method == 'wisconsin':
            data = data / data.max(axis=0)
            data = (data.T / data.sum(axis=1)).T
        else:
            normalizer = getattr(sklearn.preprocessing, method)()
            data = pd.DataFrame(
                normalizer.fit_transform(data),
                index=data.index,
                columns=data.columns
            )

        if axis == 'samples':
            data = data.T

        if not inplace:
            return data

        self.data = data

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
        main_otus = proportions.sort_values(ascending=False).index[proportions > thresh]

        return main_otus


    def select_otus(self, criteria='prevalence', n=100):
        if criteria == 'prevalence':
            scores = (self.data > 0).sum()
        elif criteria == 'abundance':
            scores = self.data.sum()
        elif criteria == 'variance':
            scores = self.data.std()
        else:
            sys.exit('Unknown criteria')

        scores.sort_values(ascending=False, inplace=True)
        return scores.index[:n]
    
    def get_most_variable_otus(self, n=100):
        proportions = self.data / self.data.sum()
        variance = proportions.sort_values(ascending=False)

        return variance.index[:n]
        

    def compute_alpha_diversity(self, metric):

        if metric == 'n_otus':
            self.alpha_diversity = self.data.sum(axis=1)
        elif metric == 'richness':
            self.alpha_diversity = (self.data > 0).sum(axis=1)
        else:
            self.alpha_diversity = self.data.apply(getattr(diversity.alpha, metric, axis=1))

        self.alpha_diversity = pd.Series(self.alpha_diversity, index=self.samples())
        self.alpha_diversity.name = metric

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
            table = handle.get('taxonomy')[:].astype(str)
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

    def clean_labels(self, inplace=True, trim=False, maxlen=20):
        data = self.data.replace(to_replace=r'[\-\(\)/]', value='_', regex=True)

        if trim:
            data = (data.fillna('')
                    .applymap(lambda x: '{}{}'.format(x[:maxlen], '...'*(len(x)>maxlen))))

        if not inplace:
            return data

        self.data = data

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

    
class SequenceData():

    def __init__(self, fasta_path=None, tree_path=None, outdir=None):
        self.fasta_path = fasta_path
        self.phylip_path = None
        self.tree_path = tree_path
        self.outdir = outdir

        if all(not x.is_file() for x in [fasta_path, tree_path] if x is not None):
            warnings.warn('No fasta or tree files found', UserWarning)

    def __repr__(self):
        return "OTUs: {} sequences".format()

    def get_otus(self):
        if self.fasta_path is not None:
            self.otus = [x.split()[0][1:] for x in open(self.fasta_path)
                         if x.startswith('>')]
        else:
            print('Not implemented yet')

    def regenerate_fasta(self, otus):
        otus = set(otus)
        filtered = [entry for entry in SeqIO.parse(self.fasta_path, 'fasta') if entry.id in otus]

        tmp_file = mkdtemp()
        self.fasta_path = Path(tmp_file, self.fasta_path.name)
        SeqIO.write(filtered, str(self.fasta_path), 'fasta')

    def update_tree(self, otus, app='FastTree', force=False):
        
        self.regenerate_fasta(otus)

        if app == 'phyml':
            phylip_path = Path(self.fasta_path.parent, "{}.phy".format(self.fasta_path.stem))
            self.tree_path = Path(self.outdir, '{}_phyml_tree.txt'.format(phylip_path.stem))
            
            if self.tree_path is not None and self.tree_path.is_file() and not force:
                print('Tree already saved in {}'.format(self.tree_path))
                return
                    
            AlignIO.convert(self.fasta_path, "fasta", phylip_path, "phylip-relaxed")

            phyml_cmd = PhymlCommandline(input=str(phylip_path))
            phyml_cmd()

        elif app == 'FastTree':
            self.tree_path = Path(self.outdir, '{}.nwk'.format(self.fasta_path.stem))

            if self.tree_path is not None and self.tree_path.is_file() and not force:
                print('Tree already saved in {}'.format(self.tree_path))
                return

            with open(str(self.tree_path), 'w') as file_handle:
                proc = subprocess.Popen(['FastTree', '-nt', str(self.fasta_path)],
                                        stdout=file_handle)
            proc.wait()

        else:
            sys.exit('Unknown application {}'.format(app))

        self.root_tree()
        
    def load_tree(self):
        return TreeNode.read(str(self.tree_path))
    
    def root_tree(self):
        tree = TreeNode.read(str(self.tree_path))
        tree = tree.root_at_midpoint()
        tree.write(str(self.tree_path))

