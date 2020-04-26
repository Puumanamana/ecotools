from itertools import product
from pathlib import Path

import pandas as pd

from bokeh.plotting import figure

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

from ecotools.api.bokeh_util import bokeh_save, bokeh_facets, bokeh_legend_out
from ecotools.api.bokeh_util import get_palette, PADDING

def pandas_to_r(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_from_pd_df = ro.conversion.py2rpy(df)

    return r_from_pd_df

def r_to_pandas(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_py = ro.conversion.rpy2py(df)

    return df_py

def run_nmds(metagenome, k=2, threads=3, trymax=200, rank='Genus', metric='braycurtis'):

    out_file = Path(metagenome.outdir, 'nmds.csv')

    if out_file.is_file():
        nmds_scores = pd.read_csv(out_file, index_col=0)
        return nmds_scores

    metagenome = metagenome.copy()

    vegan_pkg = importr('vegan')
    metagenome.group_taxa('Genus', discard_unknown=True)
    metagenome.abundance.wisconsin()
    metagenome.abundance.compute_distance_matrix(metric)
    
    dists_r = pandas_to_r(metagenome.abundance.beta_diversity.to_data_frame())

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

    metadata =  metagenome.metadata.data[metagenome.metadata.qual_vars]
    components = components.merge(metadata, left_index=True, right_index=True)

    hue_values = metadata[hue].unique()
    palette = dict(zip(hue_values, get_palette(len(hue_values))))
    
    components['color'] = [palette[x] for x in metadata[hue]]
    components.reset_index(inplace=True)

    tooltips = zip(components.drop('color', axis=1).columns, '@'+components.drop('color', axis=1).columns)

    p = figure(title="NMDS components", tooltips=list(tooltips), min_border=PADDING)

    p.circle(x='nmds_1', y='nmds_2', color='color', line_color='black',
             source=components, size=10, alpha=0.5, legend_field=hue)
    return p
