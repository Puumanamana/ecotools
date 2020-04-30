import re

import pandas as pd
import numpy as np

from ecotools.load_and_convert import calc_diversity

def subset(shared, taxonomy, metadata, alpha_div,
           sample_file=None, sample_cond=None,
           taxa_file=None, taxa_cond=None):

    sample_cond_bool = subset_samples(metadata, sample_file, sample_cond)
    taxa_cond_bool = subset_taxa(taxonomy, taxa_file, taxa_cond)
    
    shared = shared.loc[sample_cond_bool, taxa_cond_bool]
    shared = shared.loc[:, shared.sum() > 0]
    taxonomy = taxonomy.loc[shared.columns]
    metadata = metadata.loc[sample_cond_bool]

    alpha_div = alpha_div.loc[sample_cond_bool]

    if sum(~taxa_cond_bool) > 0:
        # Remove null otus
        shared = shared.loc[:, shared.sum()>0]
        taxonomy = taxonomy.loc[shared.columns]
        # We need to update the alpha diversity metrics
        alpha_div = calc_diversity(shared)

    return (shared, taxonomy, alpha_div, metadata)


# NEEDS TESTING
def subset_taxa(taxonomy, taxa_file=None, taxa_cond=None):
    """
    - taxa_file must be a csv file where the first column must be the chosen taxonomies. 
    The column title must be the rank to filter on.
    - taxa_cond is a list of conditions (only logical "and" available). Each condition is formatted as {column}@{values} where column is a taxonomic rank,  where `@` checks for inclusion in a set. {values} is a comma separated list of taxonomies.
    """

    taxa_cond_bool = np.repeat(True, taxonomy.shape[0])

    if taxa_file is not None:
        ids = pd.read_csv(taxa_file).iloc[:, 0]
        taxa_cond_bool &= np.isin(taxonomy[ids.name.title()], ids)
    if taxa_cond:
        for cond in taxa_cond:
            (col, vals) = cond.split('@')
            taxa_cond_bool &= np.isin(taxonomy[col], vals.split(','))

    return taxa_cond_bool

# NEEDS TESTING
def subset_samples(metadata, sample_file=None, sample_cond=None):
    """
    - sample_file must be a csv file without header and where the first column must be the chosen samples.
    - sample_cond is a list of conditions (only logical "and" available). Each condition is formatted as {column}{comparator}{value(s)} where column is a metadata column, comparator is in {<, >, >=, <=, ==, @}, where `@` checks for inclusion in a set. The last part is a comma separated list of values if the operator @ is chosen, and a single value otherwise.
    """

    sample_cond_bool = np.repeat(True, metadata.shape[0])
            
    if sample_file:
        ids = pd.read_csv(sample_file, header=None)[0]
        sample_cond_bool &= np.isin(metadata.index, ids)
    if sample_cond:
        for cond in sample_cond:
            (column, comparator, value) = re.findall('(.*)(@|>=|<=|<|>|==)(.*)', cond)

            column = column.strip()
            value = value.strip().replace(' ', '')

            if comparator == '@':
                sample_cond_bool &= np.isin(metadata[column], value.split(','))
            else:
                sample_cond_bool &= eval('metadata["{}"] {} {}'.format(column, comparator, value))

    return sample_cond_bool
