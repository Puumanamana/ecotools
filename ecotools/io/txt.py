from pathlib import Path

import pandas as pd
from skbio.tree import TreeNode
from Bio import SeqIO


def get_header(path, sep='\t'):
    return next(open(path)).split(sep)

def get_str_dtype(data):

    def str_len(x): return len(str(x))
    
    if len(data.shape) == 1:
        max_len = data.map(str_len).max()
    else:
        max_len = data.applymap(str_len).max().max()

    return 'S{}'.format(max_len)
    

def from_txt(path, kind=None, **kwargs):
    path = Path(path)

    if not path.is_file():
        raise FileNotFoundError(path)
    
    if kind is None:
        kind = Path(path).suffix[1:]

    parsers = {'shared': shared_from_txt,
               'taxonomy': taxonomy_from_txt,
               'metadata': metadata_from_txt}
    parsers.update({x: metadata_from_txt for x in ['csv', 'tsv']})
    parsers.update({x: tree_from_txt for x in ['nwk', 'tre']})
    parsers.update({x: fasta_from_txt for x in ['fa', 'fasta', 'fna']})

    if not kind in parsers:
        raise KeyError(f'Unknown data type {kind}')

    return parsers[kind](path, **kwargs)
       

def shared_from_txt(path, **kwargs):
    header = get_header(path)
    dtypes = dict((col, str) if col in ['Group'] else (col, int) for col in header)

    table = pd.read_csv(
        path, sep='\t', dtype=dtypes,
        keep_default_na=False, low_memory=False,
        usecols=lambda x: x not in ['label', 'numOtus']
    ).set_index('Group')
    
    table.index.name = 'SampleID'

    return table


def taxonomy_from_txt(path, species_path=None, **kwargs):

    ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

    table = (pd.read_csv(path, index_col=0, sep='\t').Taxonomy
             .str.strip(';')
             .str.split(';', expand=True))

    if species_path is not None:
        species_info = pd.read_csv(species_path, index_col=0).Species
        table = pd.concat([table, species_info], axis=1)

    table.columns = ranks[:table.shape[1]]

    return table


def metadata_from_txt(path, qual_vars=None, sep=',', **kwargs):

    if qual_vars is None:
        qual_vars = set()

    table = pd.read_csv(path, index_col=0, sep=sep,
                        dtype={col: str for col in qual_vars})
    table.index.name = 'SampleID'

    return table


def fasta_from_txt(path, **kwargs):
    data  = {x.id: x for x in SeqIO.parse(path, 'fasta')}

    return data

def tree_from_txt(path, **kwargs):
    tree = TreeNode.read(str(path))

    return tree
