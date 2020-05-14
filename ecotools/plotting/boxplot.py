import sys
import re

import numpy as np
import pandas as pd

from bokeh.plotting import figure
from bokeh.transform import dodge
from bokeh.models import HoverTool, FactorRange

from ecotools.decorators import bokeh_legend_out, bokeh_save, bokeh_facets, bokeh_cmap
from ecotools.parsing import parse_config

CFG = parse_config()['bokeh']

def boxplot(x=None, y=None, data=None, width=500, height=500, p=None, **plot_kw):

    grouped = data.groupby(x)

    quantiles = grouped[y].quantile(q=[0, 0.25, 0.5, 0.75, 1]).unstack()
    iqr = quantiles[0.75] - quantiles[0.25]

    box_data = {
        'q1': quantiles[0.25], 'q2': quantiles[0.5], 'q3': quantiles[0.75], 'iqr': iqr,
        'inf': grouped[y].agg(min), 'sup': grouped[y].agg(max),
        'upper': quantiles[0.75] + 1.5*iqr, 'lower': quantiles[0.25] - 1.5*iqr,
        'x': quantiles.index.tolist(),
    }

    if len(x) > 1:
        box_data['color'] = grouped['color'].agg('first')
        box_data[x[1]] = grouped[x[1]].agg('first')

    box_data['lower'] = pd.concat([box_data['lower'], box_data['inf']], axis=1).max(axis=1)
    box_data['upper'] = pd.concat([box_data['upper'], box_data['sup']], axis=1).min(axis=1)

    if p is None:
        p = figure(x_range=FactorRange(*box_data['x']), width=width, height=height,
                   x_axis_label=x[0], y_axis_label=y)

    leg_kw = {}
    if 'legend_field' in plot_kw:
        leg_kw = {'legend_field': plot_kw.pop('legend_field')}

    # stems
    p.segment('x', 'upper', 'x', 'q3', line_color="black", source=box_data)
    p.segment('x', 'lower', 'x', 'q1', line_color="black", source=box_data)

    # boxes
    p.vbar('x', 0.7, 'q2', 'q3', line_color='black', source=box_data, **leg_kw, **plot_kw)
    p.vbar('x', 0.7, 'q1', 'q2', line_color='black', source=box_data, **plot_kw)

    # whiskers (almost-0 height rects simpler than segments)
    p.rect('x', 'lower', 0.2, 0.01, line_color="black", source=box_data)
    p.rect('x', 'upper', 0.2, 0.01, line_color="black", source=box_data)

    p.xgrid.grid_line_color = None

    return p
    

@bokeh_save
@bokeh_facets
@bokeh_legend_out
@bokeh_cmap
def boxplot2(metagenome, x=None, hue=None, cmap=None, **kwargs):

    metagenome = metagenome.copy()

    if hue is None:
        # Make up a fake hue
        hue = ' '
        metagenome.metadata.add_var(' ', ' ')

    x_values = sorted(metagenome.metadata.factor_data(x).unique())
    hue_values = sorted(metagenome.metadata.factor_data(hue).unique())

    all_values = [(hue_i, x_i) for x_i in x_values for hue_i in hue_values]

    metric = metagenome.alpha_diversity.name
    diversity = metagenome.metadata.factor_data()
    diversity[metric] = metagenome.alpha_diversity

    if metagenome.alpha_diversity.std().max() < 1e-10:
        sys.exit('Diversity metric is constant across all levels. This usually happens with the n_otus metric if you normalized by your sample sizes. Aborting.')

    box_data = diversity.groupby([hue, x])[metric].describe().reindex(all_values)

    box_data['IQR'] = box_data['75%'] - box_data['25%']
    box_data['upper'] = np.min([box_data['max'], box_data['75%'] + 1.5*box_data['IQR']], axis=0)
    box_data['lower'] = np.max([box_data['min'], box_data['25%'] - 1.5*box_data['IQR']], axis=0)

    box_data = box_data.unstack(level=0)
    box_data.columns = [f'{col}_{hue_i}'.strip('_') for (col, hue_i) in box_data.columns]
    box_data.reset_index(inplace=True)

    scatter_data = diversity.sample(frac=1).groupby([hue, x]).head(100)

    var_tooltips = {x: re.sub("[^A-Za-z0-9_ ]", '_', x.replace("'", ""))
                    for x in metagenome.metadata.qual_vars}

    scatter_data.index.name = 'sample'
    var_tooltips['sample'] = 'sample'

    tooltips = zip(var_tooltips.values(), map(lambda x: '@{}'.format(x), var_tooltips.values()))

    scatter_data = scatter_data.rename(columns=var_tooltips)

    p = figure(x_range=x_values, min_border=CFG['padding'], tools=CFG['tools'][1:],
               plot_height=500+CFG['padding'], plot_width=1200+CFG['padding'])

    width = len(scatter_data[x].unique()) / 10
    hue_space = 1.1

    for i, hue_i in enumerate(hue_values):
        x_dodge = dodge(x, hue_space*width*(i+0.5-len(hue_values)/2), range=p.x_range)
        p.vbar(x=x_dodge, bottom=f'50%_{hue_i}', top=f'75%_{hue_i}', width=width,
               line_color="black", color=cmap[hue_i], legend_label=hue_i, source=box_data)
        p.vbar(x=x_dodge, bottom=f'25%_{hue_i}', top=f'50%_{hue_i}', width=width,
               line_color="black", color=cmap[hue_i], legend_label=hue_i, source=box_data)
        p.segment(x_dodge, f'upper_{hue_i}', x_dodge, f'75%_{hue_i}', line_color="black",
                  source=box_data)
        p.segment(x_dodge, f'lower_{hue_i}', x_dodge, f'25%_{hue_i}', line_color="black",
                  source=box_data)

        p.circle(x=x_dodge, y=metric, name='scatter', size=5,
                 line_color='black', fill_color='white', alpha=0.7,
                 source=scatter_data[scatter_data[var_tooltips[hue]]==hue_i].reset_index())

    hover = HoverTool(names=['scatter'], tooltips=list(tooltips))
    p.add_tools(hover)

    p.x_range.range_padding = 0.5
    p.x_range.factor_padding = 2 # also subgroup_padding, group_padding
    p.xaxis.major_label_orientation = 'vertical'
    p.xaxis.axis_label = x
    p.yaxis.axis_label = metric

    return p

