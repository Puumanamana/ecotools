import sys

import pandas as pd
from skbio.stats.distance import DistanceMatrix
from skbio.stats.distance import permanova, permdisp

from ecotools.decorators import pairwise, strata

@strata
@pairwise
def permanova_test(metagenome, factor, permutations=999, pairwise=False, strata=None):

    distances = DistanceMatrix(metagenome.distance_matrix)
    result = permanova(distances, factor)

    group_counts = factor.value_counts()
    group_counts.index = ['n1', 'n2']
    result = pd.concat([result, group_counts])

    return pd.DataFrame(result).T

@strata
@pairwise
def betadisper_test(metagenome, factor, permutations=999, pairwise=False, strata=None):

    distances = DistanceMatrix(metagenome.distance_matrix)
    result = permdisp(distances, factor)

    group_counts = factor.value_counts()
    group_counts.index = ['n1', 'n2']
    result = pd.concat([result, group_counts])

    return pd.DataFrame(result).T

def permutation_test(metagenome, factor, test='permanova', permutations=999,
                     pairwise=False, strata=None):

    if strata is not None and not isinstance(strata, list):
        strata = [strata]
    
    if test.lower() == 'permanova':
        results = permanova_test(metagenome, factor, permutations=permutations,
                                 pairwise=pairwise, strata=strata)
    elif test.lower() == 'betadisper':        
        results = betadisper_test(metagenome, factor, permutations=permutations,
                                  pairwise=pairwise, strata=strata)
    else:
        sys.exit('Unknown test')

    results = (results
               .rename(columns={'level_0': 'comparisons'})
               .drop(columns='level_1'))

    return results
