import os
from pathlib import Path
import re
from collections import defaultdict

from configparser import ConfigParser, ExtendedInterpolation

import h5py
import pandas as pd
import numpy as np

PARENT_DIR = Path(__file__).absolute().parent.parent

def parse_config():
    cfg = ConfigParser({'home': os.environ['HOME']},
                       interpolation=ExtendedInterpolation())
    cfg.optionxform = str
    cfg.read('{}/ikewai.ini'.format(PARENT_DIR))

    # Save parsed file for R
    with open('{}/ikewai_R.ini'.format(PARENT_DIR), 'w') as configfile:
        default_fields = dict(cfg.items('DEFAULT')).keys()
        for section in cfg.sections():
            for (name, value) in cfg.items(section):
                if name not in default_fields:
                    cfg.set(section, name, value)
        cfg.write(configfile)

    return cfg

def clean_name(s):
    s1 = re.sub('[\.\(\)\- /]', '_', s)
    s1 = re.sub('_+', '_', s1)
    s1 = re.sub('_+$', '', s1)
    s1 = s1.replace('‰', 'ratio').replace('%', 'ratio').replace('\u03b4', 'd_')

    return s1

def clean_num_vec(serie):
    if pd.api.types.is_numeric_dtype(serie):
        return serie

    serie = serie.str.replace(',', '')
    numeric = serie.apply(lambda x: str(x).replace('.', '', 1).isdigit())

    if sum(numeric) == 0:
        print("Dropping {} (non numeric column)".format(serie.name))

    serie = serie.str.replace('.*<.*', '0', regex=True)

    return pd.to_numeric(serie.where(numeric))

def load_meta(meta_file, factor_names, indices=None, include=None, exclude=None):
    meta = pd.read_csv(meta_file, index_col=0)

    for col in meta.columns:
        if col not in factor_names:
            meta[col] = clean_num_vec(meta[col])

    if include:
        for (col, value) in include.items():
            meta = meta[np.isin(meta[col], value)]
    if exclude:
        for (col, value) in exclude.items():
            meta = meta[~np.isin(meta[col], value)]

    meta.columns = [clean_name(col) for col in meta.columns]
    return meta

def load_shared(filename):
    types = defaultdict(lambda: int)
    types['Group'] = str
    shared = pd.read_csv(filename, index_col='Group', sep='\t', dtype=types,
                         keep_default_na=False, low_memory=False,
                         usecols=lambda x: x not in ['label', 'numOtus']).dropna(how='all')
    return shared

def load_tax(filename):
    tax = (pd.read_csv(filename, index_col=0, sep='\t').Taxonomy
           .str.strip(';')
           .str.split(';', expand=True))
           
    tax.columns = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus"]

    return tax

def load_diversity(filename):
    alpha_div = pd.read_csv(filename, sep='\t', index_col='group').drop('label', axis=1)

    return alpha_div

def load_data(cfg):

    shared = load_shared(cfg.get('microbial', 'shared'))
    tax = load_tax(cfg.get('microbial', 'tax'))
    alpha_div = load_diversity(cfg.get('microbial', 'alpha_div'))
    
    meta_file = cfg.get('metadata', 'file')
    factor_names = cfg.get('metadata', 'factor_names').split(',')

    include = {k: cfg.get('include', k).split(',') for k in cfg['include']
               if k not in cfg['DEFAULT']}
    exclude = {k: cfg.get('exclude', k).split(',') for k in cfg['exclude']
               if k not in cfg['DEFAULT']}
    
    meta = load_meta(meta_file, factor_names,
                     include=include,
                     exclude=exclude)

    return (shared, tax, alpha_div, meta)
        
def convert(cfg=None):
    '''
    '''

    if cfg is None:
        cfg = parse_config()

    fig_dir = Path(cfg.get('misc', 'fig_dir'))
    table_dir = Path(cfg.get('misc', 'table_dir'))

    fig_dir.mkdir(exist_ok=True, parents=True)
    table_dir.mkdir(exist_ok=True, parents=True)    
    
    (shared, tax, alpha_div, meta) = load_data(cfg)

    in_common = np.intersect1d(meta.index, shared.index)
    shared = shared.loc[in_common]
    meta = meta.loc[in_common]
    shared = shared.loc[:, shared.sum() > 0]
    alpha_div = alpha_div.loc[shared.index, :]

    h5_path = cfg.get('misc', 'h5')
    save_h5(shared, tax, alpha_div, meta, h5_path)

def save_h5(shared, tax, alpha_div, meta, output):

    samples = shared.index.values
    otu_labels = shared.columns.values

    handle = h5py.File(output, 'w')
    handle.create_dataset('abundance', data=shared.values)
    handle.create_dataset('samples', data=samples.astype('S32'))
    handle.create_dataset('OTUs', data=otu_labels.astype('S32'))
    handle.create_dataset('taxa', data=tax.loc[otu_labels].values.astype('S32'))
    handle.create_dataset('ranks', data=tax.columns.values.astype('S32'))
    handle.create_dataset('alpha_div', data=alpha_div.values)
    handle.create_dataset('alpha_metrics', data=alpha_div.columns.values.astype('S32'))    
    handle.create_dataset('metadata', data=meta.loc[samples].astype('S32').values)
    handle.create_dataset('covariates', data=meta.columns.values.astype('S32'))
    handle.close()

def load_h5(path, norm_fn=None):
    handle = h5py.File(path, 'r')

    shared = pd.DataFrame(
        handle.get('abundance')[:],
        index=handle.get('samples')[:].astype(str),
        columns=handle.get('OTUs')[:].astype(str)
    )

    taxonomy = pd.DataFrame(
        handle.get('taxa')[:].astype(str),
        index=handle.get('OTUs')[:].astype(str),
        columns=handle.get('ranks')[:].astype(str)
    ).apply(lambda x: x.str.replace('-', '_', regex=False))

    metadata = pd.DataFrame(
        handle.get('metadata')[:].astype(str),
        index=handle.get('samples')[:].astype(str),
        columns=handle.get('covariates')[:].astype(str)
    ).replace('nan', np.nan)

    alpha_div = pd.DataFrame(
        handle.get('alpha_div')[:],
        index=handle.get('samples')[:].astype(str),
        columns=handle.get('alpha_metrics')[:].astype(str)
    )

    if norm_fn is not None:
        shared = normalize(shared, fun=norm_fn)

    return (shared, taxonomy, alpha_div, metadata)

def normalize(shared, fun='assr'):
    shared = (shared.T / shared.sum(axis=1)).T

    if callable(fun):
        return fun(shared)
    elif fun == 'assr':
        return np.arcsin(np.sqrt(shared))

    print("Unknown function {}. Keeping relative abundances".format(fun))

    return shared

if __name__ == '__main__':
    convert()
