from pathlib import Path
import re

import numpy as np
import pandas as pd

from bokeh.io import output_file, save
from bokeh.plotting import figure
from bokeh.layouts import gridplot
from bokeh.models import FactorRange, Band, ColumnDataSource
from bokeh.transform import factor_cmap
from bokeh.palettes import Dark2, Accent, Category10

from ecotools.load_and_convert import parse_config
from ecotools.factors_vs_16S import parse_args

def load_results(folder, factor):
    data = []

    for f in folder.glob(f"adonis_{factor}_*.csv"):
        df = pd.read_csv(f, index_col=0)
        (lvl1, lvl2) = re.findall("_([^_]*\-vs\-[^_]*)$", f.stem)[0].split('-vs-')
        df['comparison'] = "{}-vs-{}".format(lvl1, lvl2)
        df['level_1'] = lvl1
        df['level_2'] = lvl2
        df['neg_log_padj'] = -np.log10(df['p.adjusted'])

        data.append(df)

    data = pd.concat(data)
    data.index.name = 'season'
    data = data.reset_index().set_index(['season', 'comparison'])

    full_index = pd.MultiIndex.from_product(data.index.levels)
    
    return data.reindex(full_index)

def bokeh_scatter(data, fig_dir, factor_name, threshold=0.05):
    plot_data = dict(
        x=data.index.tolist(),
        pval_adj=data['p.adjusted'].tolist(),
        neg_log10_pval_adj=data['neg_log_padj'].tolist(),
        R2=data['R2'].tolist(),
        upper=[data['neg_log_padj'].max()*1.2]*len(data),
        lower=[-np.log10(threshold)]*len(data),
        n1=data['n1'].tolist(),
        n2=data['n2'].tolist(),
    )
    plots = []

    factors = pd.MultiIndex.get_level_values(data.index, 1).unique()
    palette = Category10[10] + Dark2[8] + Accent[8]
    index_cmap = factor_cmap('x', palette=palette, factors=sorted(factors), start=1, end=2)

    tooltips = [('R2', '@R2'),
                ('adjusted p-value', '@pval_adj'),
                ('#samples in {} 1'.format(factor_name), '@n1'),
                ('#samples in {} 2'.format(factor_name), '@n2')]

    titles = {
        'neg_log10_pval_adj': 'Permanova test: -log10[adjusted p-value]',
        'R2': 'Permanova test: R2 score'
    }
    
    for metric, title in titles.items():
        p = figure(title=title, x_range=FactorRange(*plot_data['x']), tooltips=tooltips)
        p.vbar(x='x', top=metric, width=0.9, source=plot_data, line_color="white", fill_color=index_cmap)

        if metric == 'neg_log10_pval_adj':
            p.line(list(range(len(data) + len(factors))), [-np.log10(threshold)]*(len(data) + len(factors)),
                   line_color='grey', line_dash='dashed', line_width=1, line_alpha=0.5,
                   legend_label="{:.0%} significance threshold".format(threshold))

            band = Band(base='x', lower='lower', upper='upper', level='underlay', source=ColumnDataSource(plot_data),
                        line_dash='dashed', fill_alpha=0.5, line_width=1, line_color='black')
            p.add_layout(band)
            p.legend.background_fill_alpha = 0.0

        p.xaxis.major_label_orientation = "vertical"
        p.xaxis.axis_label_text_font_size = "10pt"
        plots.append(p)

    grid = gridplot(plots, ncols=1, plot_width=1500, plot_height=500)
    output_file(f"{fig_dir}/permanova_results_{factor_name}.html")
    save(grid)

if __name__ == '__main__':
    args = parse_args()
    cfg = parse_config()
    table_dir = Path(cfg.get('misc', 'table_dir'))
    fig_dir = cfg.get('misc', 'fig_dir')
    
    results = load_results(table_dir, args.factor.capitalize())

    bokeh_scatter(results, fig_dir, args.factor.lower())

