import sys
import re

import numpy as np

from bokeh.plotting import figure
from bokeh.transform import dodge
from bokeh.models import HoverTool

from ecotools.decorators import bokeh_legend_out, bokeh_save, bokeh_facets, bokeh_cmap
from ecotools.parsing import parse_config

CFG = parse_config()['bokeh']

@bokeh_save
@bokeh_facets
@bokeh_legend_out
@bokeh_cmap
def boxplot(metagenome, x=None, hue=None, cmap=None, **kwargs):

    metagenome = metagenome.copy()

    if hue is None:
        # Make up a fake hue
        hue = ' '
        metagenome.metadata.add_var(' ', ' ')

    x_values = metagenome.metadata.factor_data(x).unique()
    hue_values = metagenome.metadata.factor_data(hue).unique()

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

        p.circle(x=x_dodge, y=metric, name='scatter',
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
