from itertools import combinations, product

import pandas as pd
import numpy as np
from sklearn.metrics import pairwise_distances

from bokeh.io import save, output_file
from bokeh.palettes import inferno, Category10
from bokeh.layouts import gridplot

from load_and_convert import load_h5, parse_config
from factors_vs_chemistry import boxplot_single
from factors_vs_16S import parse_args

def calc_dist_matrix(shared, metric='braycurtis'):
    dist = pairwise_distances(shared, metric=metric, n_jobs=-1)
    dist_df = pd.DataFrame(dist, index=shared.index, columns=shared.index)
    return dist_df

def wisconsin(df):
    df = df / df.max(axis=0)
    df = df.T / df.sum(axis=1)

    return df.T

def get_normalized_distances(h5_path, factor):
    (shared, _, _, metadata) = load_h5(h5_path)

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

def boxplot_with_factor(distances, factor='Aquifer', fig_dir=None):

    distances[factor] = [g1 if g1==g2 else 'noise'
                         for (g1, g2) in distances[['g1', 'g2']].to_numpy()]

    distance_per_group = (
        distances
        .groupby(factor)
        .describe()[0]
        .reset_index()
        .rename(columns={'25%': 'q1', '50%': 'q2', '75%': 'q3', 'index': factor})
    )

    if len(distance_per_group) > 11:
        colors = inferno(n=len(distance_per_group))
    else:
        colors = Category10[len(distance_per_group)]

    cmap = {f: colors[i] for (i, f) in enumerate(distance_per_group[factor])}

    p = boxplot_single(distance_per_group, cmap, factor, 'distance')
    p.xaxis.axis_label = factor
    p.yaxis.axis_label = 'Distribution of braycurtis distances'

    output_file(f"{fig_dir}/distances_distr_by-{factor}.html")
    save(p)


def boxplot_with_factor_pairwise(distances, factor='Aquifer', fig_dir=None):
    # distances[factor] = [g1 if g1==g2 else '{}_vs_{}'.format(*sorted([g1, g2]))
    #                      for (g1, g2) in distances[['g1', 'g2']].to_numpy()]

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

if __name__ == '__main__':
    cfg = parse_config()
    args = parse_args()

    fig_dir = cfg.get('misc', 'fig_dir')

    distances = get_normalized_distances(cfg.get('misc', 'h5'), args.factor)
    boxplot_with_factor(distances, args.factor, fig_dir=fig_dir)
    boxplot_with_factor_pairwise(distances, args.factor, fig_dir=fig_dir)
