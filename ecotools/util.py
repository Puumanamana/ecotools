import re

import pandas as pd

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

def elt_or_nothing(l):
    if len(set(l)) == 1:
        return l[0]
    else:
        return None

def group_by_rank(shared, taxonomy, rank):

    rank_idx = taxonomy.columns.get_loc(rank) + 1
    ranks_to_keep = taxonomy.columns[:rank_idx]
    
    valid = ~(taxonomy[rank].str.contains('uncultured')
              | taxonomy[rank].str.contains('unclassified')
              | taxonomy[rank].str.contains('unknown'))

    shared_by_rank = shared.loc[:, valid].groupby(taxonomy[rank], axis=1).agg(sum)
    tax_by_rank = taxonomy.loc[valid].groupby(rank)[ranks_to_keep].agg(elt_or_nothing)

    if tax_by_rank[rank].isnull().sum() > 0:
        print('Some {}s have different upper level ranks'.format(rank))
        import ipdb;ipdb.set_trace()

    shared_by_rank = shared_by_rank.loc[:, shared_by_rank.sum() > 0]
    tax_by_rank = tax_by_rank.loc[shared_by_rank.columns]

    return (shared_by_rank, tax_by_rank)
    
