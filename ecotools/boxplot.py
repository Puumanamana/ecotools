import sys
import re

import numpy as np

from bokeh.plotting import figure
from bokeh.transform import dodge
from bokeh.models import HoverTool

from ecotools.bokeh_util import bokeh_save, bokeh_legend_out, bokeh_facets, get_palette
from ecotools.bokeh_util import TOOLS, PADDING

@bokeh_save
@bokeh_facets
@bokeh_legend_out
def diversity_with_meta(metagenome, columns,
                        output=None,
                        metric='richness',
                        taxa_file=None,
                        clade=False,
                        norm=True,
                        rank=None,):

    metagenome = metagenome.copy()
    metagenome.preprocess(taxa_file=taxa_file, clade=clade, norm=norm, rank=rank)
    metagenome.abundance.compute_alpha_diversity()

    x_values = metagenome.metadata.factor_data(columns[0]).unique()

    if len(columns) > 2:
        sys.exit('Too many factors. Aborting')
    
    elif len(columns) == 1:
        # Make up a fake hue
        (x_name, hue_name) = (columns[0], ' ')
        metagenome.metadata.add_var(' ', ' ')
    else:
        (x_name, hue_name) = columns

    hue_values = metagenome.metadata.factor_data(hue_name).unique()
    
    factors_combined = [(hue_i, x_i) for x_i in x_values for hue_i in hue_values]

    diversity = (metagenome.metadata.factor_data()
                 .assign(**{metric: metagenome.abundance.alpha_diversity[metric]}))

    if diversity[metric].std().max() < 1e-10:
        print('{} is constant across all {}. This usually happens with the n_otus metric if you normalized by your sample sizes.'.format(metric, 'x'.join(columns)))
        sys.exit('Aborting')

    box_data = diversity.groupby([hue_name, x_name])[metric].describe().reindex(factors_combined)

    box_data['IQR'] = box_data['75%'] - box_data['25%']
    box_data['upper'] = np.min([box_data['max'], box_data['75%'] + 1.5*box_data['IQR']], axis=0)
    box_data['lower'] = np.max([box_data['min'], box_data['25%'] - 1.5*box_data['IQR']], axis=0)

    palette = dict(zip(hue_values, get_palette(len(hue_values))))

    box_data = box_data.unstack(level=0)
    box_data.columns = [f'{col}_{hue_i}'.strip('_') for (col, hue_i) in box_data.columns]
    box_data.reset_index(inplace=True)

    scatter_data = diversity.sample(frac=1).groupby([hue_name, x_name]).head(100)

    var_tooltips = {x: re.sub("[^A-Za-z0-9_ ]", '_', x.replace("'", ""))
                    for x in metagenome.metadata.qual_vars}

    scatter_data.index.name = 'sample'
    var_tooltips['sample'] = 'sample'
        
    tooltips = zip(var_tooltips.values(), map(lambda x: '@{}'.format(x), var_tooltips.values()))

    scatter_data = scatter_data.rename(columns=var_tooltips)

    p = figure(x_range=x_values, plot_height=500+PADDING, plot_width=1200+PADDING,
               min_border=PADDING,
               tools=TOOLS[1:])

    width = len(scatter_data[x_name].unique()) / 10
    hue_space = 1.1

    for i, hue_i in enumerate(hue_values):
        x_dodge = dodge(x_name, hue_space*width*(i+0.5-len(hue_values)/2), range=p.x_range)
        p.vbar(x=x_dodge, bottom=f'50%_{hue_i}', top=f'75%_{hue_i}', width=width, line_color="black",
               source=box_data, color=palette[hue_i], legend_label=hue_i)
        p.vbar(x=x_dodge, bottom=f'25%_{hue_i}', top=f'50%_{hue_i}', width=width, line_color="black",
               source=box_data, color=palette[hue_i], legend_label=hue_i)

        p.segment(x_dodge, f'upper_{hue_i}', x_dodge, f'75%_{hue_i}', line_color="black", source=box_data)
        p.segment(x_dodge, f'lower_{hue_i}', x_dodge, f'25%_{hue_i}', line_color="black", source=box_data)

        p.circle(x=x_dodge, y=metric, name='scatter',
                 line_color='black', fill_color='white', alpha=0.7,
                 source=scatter_data[scatter_data[var_tooltips[hue_name]]==hue_i].reset_index())

    hover = HoverTool(names=['scatter'], tooltips=list(tooltips))
    p.add_tools(hover)

    p.x_range.range_padding = 0.5
    p.x_range.factor_padding = 2 # also subgroup_padding, group_padding
    p.xaxis.major_label_orientation = 'vertical'
    p.xaxis.axis_label = x_name
    p.yaxis.axis_label = metric

    return p
