import sys

import pandas as pd
from skbio.stats.distance import DistanceMatrix
from skbio.stats.distance import permanova, permdisp

from ecotools.decorators import pairwise, strata
from ecotools.rpy2_util import permanova_r

@strata
def permanova_test(metagenome, factor, permutations=999, pairwise=False, strata=None):

    factor_data = metagenome.metadata.data[factor]
    result = permanova_r(metagenome.distance_matrix, factor_data)

    group_counts = factor_data.value_counts()
    levels = (result.index
              .str.split(' vs ', expand=True)
              .to_frame()
              .applymap(lambda x: group_counts[x]))
    levels.columns = ['n1', 'n2']
    levels.index = result.index

    result = pd.concat([result, levels], axis=1)

    return result

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

    return results.reset_index()

