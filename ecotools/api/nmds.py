from itertools import product
from pathlib import Path

import pandas as pd

from bokeh.io import output_file, save
from bokeh.plotting import figure
from bokeh.layouts import gridplot
from bokeh.palettes import Category10
from bokeh.models import Legend

import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

from ecotools.group_distances import calc_dist_matrix, wisconsin
from ecotools.util import group_by_rank

def pandas_to_r(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_from_pd_df = ro.conversion.py2rpy(df)

    return r_from_pd_df

def r_to_pandas(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_py = ro.conversion.rpy2py(df)

    return df_py

def run_nmds(abundance, tax, outdir, rank='Genus', metric='braycurtis'):

    out_file = Path("{}/nmds.csv".format(outdir))

    if not out_file.is_file():

        vegan_pkg = importr('vegan')

        # df = abundance.T.groupby(tax[rank]).agg(sum).T
        (df, tax) = group_by_rank(abundance, tax, rank)
        df_norm = wisconsin(df)
        dists = calc_dist_matrix(df_norm, metric=metric)

        dists_r = pandas_to_r(dists)
        nmds_res = vegan_pkg.metaMDS(dists_r, k=3, trymax=200, parallel=60)
        nmds_scores = r_to_pandas(vegan_pkg.scores(nmds_res))
        nmds_scores = pd.DataFrame(
            nmds_scores,
            columns=["nmds_{}".format(i+1) for i in range(nmds_scores.shape[1])],
            index=abundance.index
        )

        nmds_scores.to_csv(out_file)

    else:
        nmds_scores = pd.read_csv(out_file, index_col=0)

    return nmds_scores
    

def plot_nmds(components, metadata, row=None, col=None, color=None, fig_dir=None):
 
    # tooltips = [
    #     ('Sample', '@index (@Site, @Year, @Month)'),
    #     ('Aquifer', '@Aquifer'),
    #     ('Flowpath', '@Flowpath'),
    # ]

    # info_cols = ['index', 'Site', 'Year', 'Month', 'Aquifer', 'Flowpath']
    tooltips = []

    palette = Category10[len(metadata[color].unique())]
    colormap = {st: palette[i] for i, st in enumerate(metadata[color].unique())}
    components = components.merge(metadata[color], right_index=True, left_index=True)
    components['color'] = components[color].map(lambda x: colormap[x])

    (rownames, colnames) = ([''], [''])
    if row:
        rownames = metadata[row].unique()
    if col:
        colnames = metadata[col].unique()

    (nrows, ncols) = (len(rownames), len(colnames))
        
    plots = []
    leg_space = 200

    for (i, j) in product(range(nrows), range(ncols)):
        (c1, c2) = (rownames[i], colnames[j])

        data_ij = components

        if c1:
            data_ij = data_ij[metadata[row] == c1]
        if c2:
            data_ij = data_ij[metadata[col] == c2]

        if data_ij.shape[0] < 5:
            plots.append(None)
            continue

        info_labels = ['{}: {}'.format(name, x) for name, x in zip([row, col], [c1, c2]) if x]

        p = figure(title="NMDS components ({})".format(','.join(info_labels)),
                   tooltips=tooltips,
                   min_border_right=leg_space)

        p.circle(x='nmds_1', y='nmds_2', color='color', source=data_ij, size=5, alpha=0.5,
                 legend_field=color)

        legend = Legend(items=p.legend.items.copy(), location=(10, 0))
        legend.background_fill_alpha = 0.0
        legend.border_line_alpha = 0.0

        p.legend.items.clear()
        p.add_layout(legend, 'right')
        
        plots.append(p)

    grid = gridplot(plots, ncols=ncols, plot_width=300+leg_space, plot_height=300)
    output_file("{}/nmds_by-{}-{}_color-{}.html".format(fig_dir, row, col, color))
    save(grid)
