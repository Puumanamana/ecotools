import sys

import numpy as np
import pandas as pd

from ecotools.decorators import strata
from ecotools.rpy2_util import permanova_r

@strata
def permanova_test(metagenome, factor, permutations=999, strata=None):

    factor_data = metagenome.metadata.data[factor]
    result = permanova_r(metagenome.distance_matrix, factor_data)
    result['log10_p-adj'] = -np.log10(result.pval_adj)
    result.drop(columns=['Df', 'sig', 'SumsOfSqs', 'p.value', 'pval_adj'], inplace=True)

    group_counts = factor_data.value_counts()
    levels = (result.index
              .str.split(' vs ', expand=True)
              .to_frame()
              .applymap(lambda x: group_counts[x]))
    levels.columns = ['n1', 'n2']
    levels.index = result.index

    result = pd.concat([result, levels], axis=1)

    return result

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

