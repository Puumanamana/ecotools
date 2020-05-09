from itertools import combinations

import pandas as pd
from bokeh.plotting import figure

from ecotools.parsing import parse_config
# from ecotools.util import get_palette
from ecotools.decorators import bokeh_save, bokeh_facets, bokeh_legend_out, bokeh_cmap

CFG = parse_config()['bokeh']

@bokeh_save
@bokeh_facets
@bokeh_legend_out
@bokeh_cmap
def scatterplot(metagenome, values, hue=None, size=10,
                width=400, height=400, labels=None, **kwargs):

    if labels is None:
        labels = [x for x in values.columns if x not in kwargs.values()]
    
    metadata =  metagenome.metadata.data[metagenome.metadata.qual_vars]
    data = values[labels].merge(metadata, left_index=True, right_index=True).sort_values(by=hue)

    data['color'] = [kwargs['cmap'][x] for x in data[hue]]
    data.reset_index(inplace=True)

    tooltips = zip(data.drop('color', axis=1).columns, '@'+data.drop('color', axis=1).columns)

    plots = []
    for (label1, label2) in combinations(labels, 2):

        p = figure(x_axis_label=label1, y_axis_label=label2,
                   width=width, height=height,
                   tooltips=list(tooltips), min_border=CFG['padding'])
        p.circle(x=label1, y=label2, color='color', line_color='gray',
                 source=data, size=size, alpha=0.5, legend_field=hue)
        plots.append(p)

    return plots

@bokeh_save
@bokeh_cmap
def corrplot(metagenome, hue=None, covariates=None, other=None,
             width=400, height=400, size=10, **kwargs):

    data = []
    if covariates is not None:
        data.append(metagenome.metadata.data[covariates])
    if other is not None:
        data.append(other)

    data = pd.concat(data, axis=1).select_dtypes(include='number')

    if hue is not None:
        data['color'] = [kwargs['cmap'][x] for x in metagenome.metadata.factor_data(hue)]
    
    data.reset_index(inplace=True)

    tooltips = zip(data.drop('color', axis=1).columns, '@'+data.drop('color', axis=1).columns)

    plots = []
    for (label1, label2) in combinations(data.columns, 2):
        p = figure(x_axis_label=label1, y_axis_label=label2,
                   width=width, height=height,
                   tooltips=list(tooltips), min_border=CFG['padding'])
        p.circle(x=label1, y=label2, color='color', line_color='gray',
                 source=data, size=size, alpha=0.5, legend_field=hue)
        plots.append(p)

    return plots




    
