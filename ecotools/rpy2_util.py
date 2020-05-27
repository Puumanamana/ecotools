from itertools import combinations
import pandas as pd

import rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

rpy2.rinterface_lib.callbacks.consolewrite_warnerror = lambda x: None

def pandas_to_r(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_from_pd_df = ro.conversion.py2rpy(df)

    return r_from_pd_df

def r_to_pandas(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_py = ro.conversion.rpy2py(df)

    return df_py

def pandas_to_distR(dist_df):
    stats_pkg = importr('stats')
    dists_r = stats_pkg.dist(pandas_to_r(dist_df))

    return dists_r
    

def vegdist(abundance, metric='bray', binary=False, r_obj=False):
    if metric.lower().startswith('bray'):
        metric = 'bray'
    vegan_pkg = importr('vegan')
    abundance_r = pandas_to_r(abundance.sort_index())
    
    dists_r = vegan_pkg.vegdist(abundance_r, method=metric,
                                binary=binary, diag=False)

    if r_obj:
        return dists_r

    indices = combinations(abundance.index, 2)
    dists = pd.Series(r_to_pandas(dists_r), index=pd.Index(indices))
    diag = pd.Series(1, index=pd.Index(zip(abundance.index, abundance.index)))
    
    dists = pd.concat([dists, diag])
    return dists.sort_index()

def metamds(dists_r, metric='bray', k=2, trymax=500, parallel=20):
    vegan_pkg = importr('vegan')
    nmds_obj = vegan_pkg.metaMDS(dists_r, distance=metric, k=k,
                                 trymax=trymax, parallel=parallel)
    nmds_components = vegan_pkg.scores(nmds_obj)

    nmds_components = pd.DataFrame(
        r_to_pandas(nmds_components),
        columns=["nmds_{}".format(i+1) for i in range(k)],
    )

    return nmds_components

def pcoa_ape(dists_r, k=2):

    ape_pkg = importr('ape')
    pcoa_obj = ape_pkg.pcoa(dists_r)

    components = pd.DataFrame(
        r_to_pandas(pcoa_obj.rx2('vectors'))[:, :k],
        columns=["pcoa_{}".format(i+1) for i in range(k)],
    )

    return components

def permanova_r(dists_r, factors, permutations=9999):
    factors_r = pandas_to_r(factors)
     
    pwAdonis_pkg = importr('pairwiseAdonis')

    r_model = pwAdonis_pkg.pairwise_adonis(dists_r, factors_r,
                                           perm=permutations)

    r_res_items = list(r_model.items())
    index = pd.Index(r_res_items[0][1].levels, name=r_res_items[0][0])

    model_results = pd.DataFrame(
        dict(r_res_items[1:]),
        columns=r_model.colnames[1:],
        index=index
    ).rename(columns={'p.adjusted': 'pval_adj', 'F.Model': 'statistic'})     

    return model_results
