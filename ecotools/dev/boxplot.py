import sys
import numpy as np
import pandas as pd

from bokeh.plotting import figure
from bokeh.transform import jitter
from bokeh.models import HoverTool

from util import bokeh_save, get_palette, TOOLS

@bokeh_save
def diversity_with_meta(metagenome, column, output=None,
                        metric='richness',
                        taxa_file=None,
                        norm=True,
                        rank=None):

    metagenome = metagenome.copy()
    metagenome.preprocess(taxa_file=taxa_file, norm=norm, rank=rank)
    metagenome.abundance.calc_diversity()

    diversity = (metagenome.metadata.factor_data()
                 .assign(**{metric: metagenome.abundance.alpha_diversity[metric]}))

    if diversity[metric].std().max() < 1e-10:
        print('{} is constant across all {}. This usually happens with the n_otus metric if you normalized by your sample sizes.'.format(metric, column))
        sys.exit('Aborting')

    box_data = diversity.groupby(column)[metric].describe()

    box_data['IQR'] = box_data['75%'] - box_data['25%']
    box_data['upper'] = np.min([box_data['max'], box_data['75%'] + 1.5*box_data['IQR']], axis=0)
    box_data['lower'] = np.max([box_data['min'], box_data['25%'] - 1.5*box_data['IQR']], axis=0)

    factors = box_data.index.tolist()
    if all(x.isdigit() for x in factors):
        factors = sorted(factors, key=int)

    palette = dict(zip(factors, get_palette(len(factors))))

    box_data['color'] = [palette[g] for g in factors]
    box_data.reset_index(inplace=True)

    p = figure(title=output, tools=TOOLS[1:],
               x_range=factors)

    # stems
    p.segment(column, 'upper', column, '75%', line_color="black", source=box_data)
    p.segment(column, 'lower', column, '25%', line_color="black", source=box_data)

    # boxes
    p.vbar(x=column, width=0.7, bottom='50%', top='75%',
           line_color="black", color='color', source=box_data)
    p.vbar(x=column, width=0.7, bottom='25%', top='50%',
           line_color="black", color='color', source=box_data)

    # scatter points
    scatter_data = diversity.sample(frac=1).groupby(column).head(500)

    scatter = p.circle(x=jitter(column, 0.2, range=p.x_range), y=metric,
                       line_color='black', fill_color='white', alpha=0.8,
                       source=scatter_data)

    tooltips = zip(metagenome.metadata.qual_vars, '@'+metagenome.metadata.qual_vars)
    hover = HoverTool(renderers=[scatter], tooltips=list(tooltips))
    p.add_tools(hover)

    p.xaxis.major_label_orientation = 'vertical'
    p.xaxis.axis_label = column
    p.yaxis.axis_label = metric

    return p
