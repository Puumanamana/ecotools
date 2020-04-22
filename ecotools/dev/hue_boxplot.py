import sys
import numpy as np

from bokeh.plotting import figure
from bokeh.transform import dodge
from bokeh.models import HoverTool

from util import bokeh_save, bokeh_legend_out, get_palette, TOOLS

@bokeh_save
@bokeh_legend_out
def diversity_with_meta(metagenome, columns,
                        output=None,
                        metric='richness',
                        taxa_file=None,
                        norm=True,
                        rank=None,
                        padding=200):

    metagenome = metagenome.copy()
    metagenome.preprocess(taxa_file=taxa_file, norm=norm, rank=rank)
    metagenome.abundance.calc_diversity()

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

    # Properly sort factors levels
    if all(x.isdigit() for x in x_values):
        x_values = sorted(x_values, key=float)
    if all(x.isdigit() for x in hue_values):
        hue_values = sorted(hue_values, key=float)
    
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

    scatter_data = diversity.sample(frac=1).groupby([hue_name, x_name]).head(500)

    p = figure(x_range=x_values, plot_height=500+padding, plot_width=1200+padding,
               min_border=padding,
               title=output.name, tools=TOOLS[1:])

    width = 0.8 / len(hue_values)

    for i, hue_i in enumerate(hue_values):
        x_dodge = dodge(x_name, 1.1*width*(i+0.5-len(hue_values)/2), range=p.x_range)
        p.vbar(x=x_dodge, bottom=f'50%_{hue_i}', top=f'75%_{hue_i}', width=width, line_color="black",
               source=box_data, color=palette[hue_i], legend_label=hue_i)
        p.vbar(x=x_dodge, bottom=f'25%_{hue_i}', top=f'50%_{hue_i}', width=width, line_color="black",
               source=box_data, color=palette[hue_i], legend_label=hue_i)

        p.segment(x_dodge, f'upper_{hue_i}', x_dodge, f'75%_{hue_i}', line_color="black", source=box_data)
        p.segment(x_dodge, f'lower_{hue_i}', x_dodge, f'25%_{hue_i}', line_color="black", source=box_data)

        p.circle(x=x_dodge, y=metric, name='scatter',
                 line_color='black', fill_color='white', alpha=0.7,
                 source=scatter_data[scatter_data[hue_name]==hue_i])

    tooltips = zip(metagenome.metadata.qual_vars, '@'+metagenome.metadata.qual_vars)
    hover = HoverTool(names=['scatter'], tooltips=list(tooltips))
    p.add_tools(hover)

    p.xaxis.major_label_orientation = 'vertical'
    p.xaxis.axis_label = x_name
    p.yaxis.axis_label = metric

    return p
