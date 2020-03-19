from itertools import product

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from bokeh.io import output_file, save
from bokeh.plotting import figure
from bokeh.layouts import gridplot
from bokeh.palettes import Category10

from load_and_convert import load_h5, parse_config

TOOLS = ['hover', 'box_zoom']

DIV_ID = ['sobs', 'shannoneven']
# COVARS = [clean_name(x) for x in ['pH', 'Potassium', 'Sulfate', 'Sodium', 'Magnesium', 'Chloride', 'Calcium', 'Temp (C)', 'DO (mg/L)', 'Si', 'Sr', 'Sp. Cond. (uS/cm)']]
# VARS = DIV_ID + COVARS

def pairwise_scatter(shared, tax, alphaD, meta, factor='Aquifer', fig_dir=None):

    # meta[COVARS] = meta[COVARS].astype(float)
    chem_cols = meta.columns[meta.dtypes=='float']
    chem_cols = np.setdiff1d(chem_cols, DIV_ID)

    meta[chem_cols] = (meta[chem_cols] - meta[chem_cols].mean()) / meta[chem_cols].std()

    data = shared.T.groupby(tax.Genus).agg(sum).T
    data = pd.concat([data, meta, alphaD.astype(float)], axis=1)
    data = data[~meta[factor].isnull()].reset_index()

    covars = chem_cols + DIV_ID

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

        corr = '{:.2f}'.format(pearsonr(df[c1], df[c2])[0])

        p = figure(title=f"{c1} vs {c2} (corr={corr})", tooltips=tooltips, tools=TOOLS,
                   min_border_right=leg_space)
        p.circle(x=c1, y=c2, color='color', source=df, size=5, alpha=0.5, hover_color="red")

        plots.append(p)

    grid = gridplot(plots, ncols=len(covars), plot_width=300+leg_space, plot_height=300)
    output_file(f"{fig_dir}/covariates_correlation_by-{factor}.html")
    save(grid)

if __name__ == '__main__':
    cfg = parse_config()
    
    (shared, taxonomy, alpha_div, metadata) = load_h5(cfg.get('misc', 'h5'))
    
    print('Data loaded')

    pairwise_scatter(shared, taxonomy, alpha_div, metadata, factor='Flowpath', fig_dir=cfg.get('misc', 'fig_dir'))
    pairwise_scatter(shared, taxonomy, alpha_div, metadata, factor='Aquifer', fig_dir=cfg.get('misc', 'fig_dir'))
