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
    nmds_obj = vegan_pkg.metaMDS(dists_r, k=k, trymax=trymax, parallel=parallel)
    nmds_components = vegan_pkg.scores(nmds_obj)

    nmds_components = pd.DataFrame(
        r_to_pandas(nmds_components),
        columns=["nmds_{}".format(i+1) for i in range(k)],
        index=dists.index
    )

    return nmds_components

