import pandas as pd

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

def pandas_to_r(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_from_pd_df = ro.conversion.py2rpy(df)

    return r_from_pd_df

def r_to_pandas(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_py = ro.conversion.rpy2py(df)

    return df_py

def metamds(dists, k=3, trymax=200, parallel=3):
    vegan_pkg = importr('vegan')
    dists_r = pandas_to_r(dists)
    nmds_obj = vegan_pkg.metaMDS(dists_r, k=k, trymax=trymax, parallel=parallel, trace=False)
    nmds_components = vegan_pkg.scores(nmds_obj)

    nmds_components = pd.DataFrame(
        r_to_pandas(nmds_components),
        columns=["nmds_{}".format(i+1) for i in range(k)],
        index=dists.index
    )

    return nmds_components


def permanova_r(distances, factors, permutations=9999):
     pwAdonis_pkg = importr('pairwiseAdonis')
     dists_r = pandas_to_r(distances)
     factors_r = pandas_to_r(factors)

     r_model = pwAdonis_pkg.pairwise_adonis(dists_r, factors_r,
                                            perm=permutations)

     r_res_items = list(r_model.items())
     index = pd.Index(r_res_items[0][1].levels, name = r_res_items[0][0])

     model_results = pd.DataFrame(
         dict(r_res_items[1:]),
         columns=r_model.colnames[1:],
         index=index
     ).rename(columns={'p.adjusted': 'pval_adj'})

     return model_results
