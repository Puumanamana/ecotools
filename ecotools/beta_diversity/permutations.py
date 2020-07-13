import sys

import numpy as np
import pandas as pd

from ecotools.decorators import strata, timer
from ecotools.rpy2_util import permanova_r, pandas_to_distR

@strata
def permanova_test(metagenome, factor, permutations=999, strata=None):

    metagenome.compute_distance_matrix(cache=True, strata=strata)
    factor_data = metagenome.metadata.data[factor]
    
    dists_r = pandas_to_distR(metagenome.distance_matrix.unstack())
    result = permanova_r(dists_r, factor_data)
    result['neg_log10_padj'] = -np.log10(result.pval_adj)
    result.drop(columns=['Df', 'sig', 'SumsOfSqs', 'p.value'], inplace=True)

    group_counts = factor_data.value_counts()
    levels = (result.index
              .str.split(' vs ', expand=True)
              .to_frame()
              .applymap(lambda x: group_counts[x]))
    levels.columns = ['n1', 'n2']
    levels.index = result.index

    result = pd.concat([result, levels], axis=1)

    return result

@timer
def permutation_test(metagenome, factor, test='permanova', permutations=999,
                     strata=None):

    if strata is not None and not isinstance(strata, list):
        strata = [strata]
    
    if test.lower() == 'permanova':
        results = permanova_test(metagenome, factor, permutations=permutations,
                                 strata=strata)
    else:
        sys.exit('Not implemented')

    return results.reset_index()

