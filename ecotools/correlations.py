from itertools import product

import pandas as pd
from scipy.stats import pearsonr

from bokeh.io import output_file, save
from bokeh.plotting import figure
from bokeh.layouts import gridplot
from bokeh.palettes import Category10

from ecotools.util import group_by_rank

TOOLS = ['hover', 'box_zoom']
DIV_ID = ['sobs', 'shannoneven']

def pairwise_scatter(shared, tax, alphaD, meta, factor='Aquifer', fig_dir=None):

    meta_quant = meta.select_dtypes(include='number')
    quant_counts = meta_quant.count().sort_values(ascending=False)
    
    meta[meta_quant.columns] = (meta_quant - meta_quant.mean()) / meta_quant.std()

    (data, tax) = group_by_rank(shared, tax, 'Genus')
    data = pd.concat([data, meta, alphaD.astype(float)], axis=1)
    data = data[~meta[factor].isnull()].reset_index()

    covars = DIV_ID + quant_counts.index.tolist()[:20]

    tooltips = [
        ('Sample', '@index (@Site, @Year, @Month)'),
        ('Aquifer', '@Aquifer'),
        ('Flowpath', '@Flowpath'),
    ]

    info_cols = ['index', 'Site', 'Year', 'Month', 'Aquifer', 'Flowpath']
    
    palette = Category10[len(data[factor].unique())]
    colormap = {st: palette[i] for i, st in enumerate(data[factor].unique())}
    
    plots = []
    leg_space = 100

    for i, j in product(range(len(covars)), repeat=2):

        (c1, c2) = (covars[i], covars[j])
        
        if i == 1 and j == 1:
            # Make legend on an empty space
            p = figure(outline_line_color=None)
            p.grid.visible = False
            p.axis.visible = False
            
            [p.circle(x=[0], y=[0], color=color, legend_label=label)
             for label, color in colormap.items()]
            
            plots.append(p)
            continue

        if i >= j:
            plots.append(None)
            continue

        df = data.loc[~data[[c1, c2]].duplicated(), info_cols+[c1, c2]].dropna(how='any')
        df.loc[:, 'color'] = [colormap[st] for i, st in enumerate(df[factor])]

        if df.shape[0] < 5:
            plots.append(None)
            continue

        corr = '{:.2f}'.format(pearsonr(df[c1], df[c2])[0])

        p = figure(title=f"{c1} vs {c2} (corr={corr})", tooltips=tooltips, tools=TOOLS,
                   min_border_right=leg_space)
        p.circle(x=c1, y=c2, color='color', source=df, size=5, alpha=0.5, hover_color="red")

        plots.append(p)

    grid = gridplot(plots, ncols=len(covars), plot_width=300+leg_space, plot_height=300)
    output_file(f"{fig_dir}/covariates_correlation_by-{factor}.html")
    save(grid)
