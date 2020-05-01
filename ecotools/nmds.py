from itertools import combinations
from pathlib import Path

import pandas as pd

from bokeh.plotting import figure

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

from ecotools.bokeh_util import bokeh_save, bokeh_facets, bokeh_legend_out
from ecotools.bokeh_util import get_palette, PADDING

def pandas_to_r(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_from_pd_df = ro.conversion.py2rpy(df)

    return r_from_pd_df

def r_to_pandas(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_py = ro.conversion.rpy2py(df)

    return df_py

def nmds(metagenome, k=2, threads=3, trymax=200, cache=False,
         rank='Genus', metric='braycurtis'):

    out_file = Path(metagenome.outdir, 'nmds_{}_on_{}_k={}.csv'.format(metric, rank.lower(), k))

    if out_file.is_file() and cache:
        nmds_scores = pd.read_csv(out_file, index_col=0)
        return nmds_scores

    metagenome = metagenome.copy()

    vegan_pkg = importr('vegan')
    
    dists_r = pandas_to_r(metagenome.abundance.beta_diversity)

    nmds_res = vegan_pkg.metaMDS(dists_r, k=k, trymax=trymax, parallel=threads)
    nmds_scores = r_to_pandas(vegan_pkg.scores(nmds_res))
    nmds_scores = pd.DataFrame(
        nmds_scores,
        columns=["nmds_{}".format(i+1) for i in range(nmds_scores.shape[1])],
        index=metagenome.samples()
    )

    nmds_scores.to_csv(out_file)

    return nmds_scores    

@bokeh_save
@bokeh_facets
@bokeh_legend_out
def plot_nmds(metagenome, components, hue=None, output=None):

    n_components = components.shape[1]
    metadata =  metagenome.metadata.data[metagenome.metadata.qual_vars]
    components = components.merge(metadata, left_index=True, right_index=True)

    hue_values = metadata[hue].unique()
    palette = dict(zip(hue_values, get_palette(len(hue_values))))
    
    components['color'] = [palette[x] for x in metadata[hue]]
    components.reset_index(inplace=True)

    tooltips = zip(components.drop('color', axis=1).columns, '@'+components.drop('color', axis=1).columns)

    plots = []
    for (k1, k2) in combinations(range(1, n_components+1), 2):
        (label1, label2) = ('nmds_{}'.format(k1), 'nmds_{}'.format(k2))

        p = figure(title="NMDS components ({} vs {})".format(k1, k2),
                   tooltips=list(tooltips), min_border=PADDING)
        p.circle(x=label1, y=label2, color='color', line_color='gray',
                 source=components, size=10, alpha=0.5, legend_field=hue)
        plots.append(p)

    return plots
