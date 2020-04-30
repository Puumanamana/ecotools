from itertools import combinations, product

import pandas as pd
import numpy as np
from sklearn.metrics import pairwise_distances

from bokeh.plotting import figure
from bokeh.io import save, output_file
from bokeh.palettes import inferno, Category10
from bokeh.layouts import gridplot
from bokeh.transform import jitter

def calc_dist_matrix(shared, metric='braycurtis'):
    dist = pairwise_distances(shared, metric=metric, n_jobs=-1)
    dist_df = pd.DataFrame(dist, index=shared.index, columns=shared.index)
    return dist_df

def wisconsin(df):
    df = df / df.max(axis=0)
    df = df.T / df.sum(axis=1)

    return df.T

def get_normalized_distances(shared, metadata, factor):
    print('Normalization')
    shared_norm = wisconsin(np.sqrt(shared))
    print('Distance matrix')
    distances = calc_dist_matrix(shared_norm).stack()
    combs = list(combinations(shared.index, 2))

    distances = distances.loc[combs].reset_index()

    # Annotate comparisons
    factors = metadata[factor]
    distances['g1'] = factors.loc[distances.level_0].to_numpy()
    distances['g2'] = factors.loc[distances.level_1].to_numpy()

    distances = distances[~distances.g1.isnull() & ~distances.g2.isnull()]

    return distances

def boxplot_with_factor_pairwise(distances, factor='Aquifer', fig_dir=None):

    plots = []
    for l1, l2 in product(sorted(set(distances.g1)), repeat=2):
        if l1 >= l2:
            plots.append(None)
            continue

        comparisons = []

        for (g1, g2, dist) in distances[['g1', 'g2', 0]].to_numpy():
            if {g1, g2} == {l1, l2} or {g1, g2} == {l1} or {g1, g2} == {l2}:
                group = {g1, g2}.pop() if g1 == g2 else f'{l1}_vs_{l2}'
                comparisons.append([l1, l2, dist, group])

        comparisons = pd.DataFrame(comparisons, columns=['g1', 'g2', 'dist', factor])

        df = (
            comparisons
            .groupby(factor)
            .describe()['dist']
            .reset_index()
            .rename(columns={'25%': 'q1', '50%': 'q2', '75%': 'q3', 'index': factor})
        )

        if len(df) > 11:
            colors = inferno(n=len(df))
        else:
            colors = Category10[len(df)]

        cmap = {f: colors[i] for (i, f) in enumerate(df[factor])}

        scatter_data = comparisons.drop(['g1', 'g2'], axis=1).rename(columns={'dist': 'y'})

        p = boxplot_single(df, cmap, factor, 'distance', scatter_data=scatter_data)

        plots.append(p)

    grid = gridplot(plots, ncols=len(set(distances.g1)))

    output_file(f"{fig_dir}/pairwise_distances_distr_by-{factor}.html")
    save(grid)

    return


def boxplot_single(info, cmap, col1, col2=None, scatter_data=None):

    info['IQR'] = info.q3 - info.q1

    info['upper'] = np.min([info['max'], info.q3 + 1.5*info.IQR], axis=0)
    info['lower'] = np.max([info['min'], info.q1 - 1.5*info.IQR], axis=0)

    info['color'] = [cmap[f] for f in info[col1]]

    tooltips = []
    if scatter_data is not None:
        tooltips = scatter_data.drop(['y', col1], axis=1).columns
        tooltips = zip(tooltips, '@'+tooltips)
    
    p = figure(title=f"{col2} vs {col1}",
               background_fill_color="#efefef",
               plot_width=300, plot_height=400,
               tooltips=list(tooltips),
               x_range=info[col1].tolist())

    # stems
    p.segment(col1, 'upper', col1, 'q3', line_color="black", source=info)
    p.segment(col1, 'lower', col1, 'q1', line_color="black", source=info)

    # boxes
    p.vbar(x=col1, width=0.7, bottom='q2', top='q3',
           line_color="black", color='color', source=info)
    p.vbar(x=col1, width=0.7, bottom='q1', top='q2',
           line_color="black", color='color', source=info)

    if scatter_data is not None:
        scatter_data = scatter_data.sample(frac=1).groupby(col1).head(500)
        scatter_data['color'] = [cmap[x] for x in scatter_data[col1]]

        p.circle(x=jitter(col1, 0.2, range=p.x_range), y='y',
                 line_color='black', fill_color='color', alpha=0.5,
                 source=scatter_data)

    # # whiskers (almost-0 height rects simpler than segments)
    # h = np.abs(info.q3).mean()
    # p.rect(x=col1, y='lower', width=0.2, height=0.01*h, color="black", source=info)
    # p.rect(x=col1, y='upper', width=0.2, height=0.01*h, color="black", source=info)

    p.xaxis.major_label_orientation = "vertical"

    return p
