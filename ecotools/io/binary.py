import sys
from pathlib import Path

import h5py

import pandas as pd

from ecotools.core.abundance_table import AbundanceTable
from ecotools.core.taxonomy_table import TaxonomyTable
from ecotools.core.metadata_table import MetadataTable

def from_h5(path, kind=None):
    path = Path(path)
    
    if not path.is_file():
        raise FileNotFoundError(path)
    
    parsers = {'abundance': shared_from_h5,
               'taxonomy': taxonomy_from_h5,
               'metadata': metadata_from_h5}

    if not kind in parsers:
        raise KeyError(f'Unknown data type {kind}')

    return parsers[kind](path)
       

def shared_from_h5(path):
    with h5py.File(path, 'r') as handle:
        table = handle.get('abundance')[:]
        samples = handle.get('samples')[:].astype(str)
        OTUs = handle.get('OTUs')[:].astype(str)

    table = AbundanceTable(
        data=pd.DataFrame(table, index=samples, columns=OTUs),
    )
    
    return table


def taxonomy_from_h5(path):
    with h5py.File(path, 'r') as handle:
        table = handle.get('taxonomy')[:].astype(str)
        levels = handle.get('tax_levels')[:].astype(str)
        OTUs = handle.get('OTUs')[:].astype(str)

    table = TaxonomyTable(
        data=pd.DataFrame(table, index=OTUs, columns=levels)
    )

    return table

def metadata_from_h5(path):
    with h5py.File(path, 'r') as handle:
        table_qual = handle.get('metadata_qual')[:].astype(str)
        table_quant = handle.get('metadata_quant')[:]
        samples = handle.get('samples')[:].astype(str)
        
        quant_vars = handle.get('quant_vars')[:].astype(str)
        qual_vars = handle.get('qual_vars')[:].astype(str)

    table = pd.concat([
        pd.DataFrame(table_qual, index=samples, columns=qual_vars),
        pd.DataFrame(table_quant, index=samples, columns=quant_vars),
    ], axis=1)

    table = MetadataTable(
        data=table,
    )
    
    return table


def fasta_from_h5(path):
    sys.exit('Fasta from h5 not implentemented yet')

def tree_from_h5(path):
    sys.exit('Tree from h5 not implentemented yet')
