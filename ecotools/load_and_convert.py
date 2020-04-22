from pathlib import Path
from collections import defaultdict

import h5py
import pandas as pd
import numpy as np
from skbio import diversity

from ecotools.parser import parse_config, get_strata
from ecotools.util import clean_name, clean_num_vec

def load_meta(meta_file, factor_names, indices=None, exclude=None, min_samples=0):

    meta = pd.read_csv(meta_file, index_col=0)

    for col in meta.columns:
        if col in factor_names:
            meta[col] = meta[col].astype(str).str.replace('.0$', '', regex=False)
        else:
            meta[col] = clean_num_vec(meta[col])

    if exclude:
        for (col, value) in exclude.items():
            meta = meta[~np.isin(meta[col], value)]

    meta.dropna(thresh=min_samples, inplace=True, axis=1)

    # remove numeric columns with std = 0
    quant_std = meta.select_dtypes(include='number').std()
    meta.drop(quant_std.index[quant_std == 0], axis=1, inplace=True)

    meta.columns = [clean_name(col) for col in meta.columns]
    import ipdb;ipdb.set_trace()

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

def calc_diversity(shared):
    alpha_div = pd.DataFrame({
        'n_OTUs': shared.sum(axis=1),
        'richness': (shared > 0).sum(axis=1),
        'chao': shared.apply(diversity.alpha.chao1, axis=1),
        # 'chao_ci': shared.apply(diversity.alpha.chao1_ci, axis=1),
        'shannon': shared.apply(diversity.alpha.shannon, axis=1),
        # 'shannon_even': ,
    })

    return alpha_div

def load_data(cfg):

    shared = load_shared(cfg.get('microbial', 'shared'))
    tax = load_tax(cfg.get('microbial', 'tax'))
    alpha_div = calc_diversity(shared)
    
    meta_file = cfg.get('metadata', 'file')
    factor_names = cfg.get('metadata', 'factor_names').split(',')
    min_samples = int(cfg.get('metadata', 'min_samples'))

    exclude = {k: cfg.get('exclude', k).split(',') for k in cfg['exclude']
               if k not in cfg['DEFAULT']}
    
    meta = load_meta(meta_file, factor_names,
                     exclude=exclude,
                     min_samples=min_samples)

    return (shared, tax, alpha_div, meta)
        
def convert(cfg=None):
    '''
    '''

    if cfg is None:
        cfg = parse_config()

    (shared, tax, alphadiv, meta) = load_data(cfg)

    for comb_name, comb in get_strata(cfg).items():
        if not comb_name:
            conditions = np.repeat(True, len(meta))
        else:
            conditions = np.prod([meta[col].astype(str) == val for (col, val) in comb],
                                 axis=0).astype(bool)

        in_common = np.intersect1d(meta[conditions].index, shared.index)

        shared_sub = shared.loc[in_common]
        meta_sub = meta.loc[in_common]
        shared_sub = shared_sub.loc[:, shared_sub.sum() > 0]
        alphadiv_sub = alphadiv.loc[in_common, :]
        tax_sub = tax.loc[shared_sub.columns]

        h5_path = Path(cfg.get('misc', 'h5'))

        save_h5(shared_sub, tax_sub, alphadiv_sub, meta_sub,
                output="{}_{}.h5".format(h5_path, comb_name))

def save_h5(shared, tax, alpha_div, meta, output=None):

    meta_quant = meta.select_dtypes(include='number')
    meta_qual = meta.drop(meta_quant.columns, axis=1)

    samples = shared.index.to_numpy()
    otu_labels = shared.columns.to_numpy()

    handle = h5py.File(output, 'w')
    handle.create_dataset('abundance', data=shared.to_numpy())
    handle.create_dataset('samples', data=samples.astype('S32'))
    handle.create_dataset('OTUs', data=otu_labels.astype('S32'))
    handle.create_dataset('taxa', data=tax.loc[otu_labels].to_numpy().astype('S32'))
    handle.create_dataset('ranks', data=tax.columns.to_numpy().astype('S32'))
    handle.create_dataset('alpha_div', data=alpha_div.to_numpy())
    handle.create_dataset('alpha_metrics', data=alpha_div.columns.to_numpy().astype('S32'))

    handle.create_dataset('metadata_quant', data=meta_quant.loc[samples].to_numpy())
    handle.create_dataset('metadata_qual', data=meta_qual.loc[samples].astype('S32').to_numpy())
    handle.create_dataset('quant_vars', data=meta_quant.columns.to_numpy().astype('S32'))
    handle.create_dataset('qual_vars', data=meta_qual.columns.to_numpy().astype('S32'))    
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

    meta_quant = pd.DataFrame(
        handle.get('metadata_quant')[:],
        index=handle.get('samples')[:].astype(str),
        columns=handle.get('quant_vars')[:].astype(str)
    )

    meta_qual = pd.DataFrame(
        handle.get('metadata_qual')[:].astype(str),
        index=handle.get('samples')[:].astype(str),
        columns=handle.get('qual_vars')[:].astype(str)
    ).replace('nan', np.nan)

    metadata = pd.concat([meta_qual, meta_quant], axis=1)

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
