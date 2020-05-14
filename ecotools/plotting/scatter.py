from itertools import combinations

import pandas as pd
from bokeh.plotting import figure
from bokeh.transform import jitter
from bokeh.models import FactorRange

from ecotools.parsing import parse_config
from ecotools.plotting.facetgrid import BokehFacetGrid

CFG = parse_config()['bokeh']

def scatter(x=None, y=None, hue=None, data=None, width=500, height=500, p=None, **plot_kw):

    if hue is not None:
        plot_kw['legend_field'] = hue
        plot_kw['fill_color'] = 'color'

    if p is None:
        p = figure(width=width, height=height, min_border=100)

    p.circle(x=x, y=y, line_color='gray', fill_color='color',
             x_axis_label=x, y_axis_label=y,
             source=data, size=10, alpha=0.5, **plot_kw)

    return p

def swarmplot(x=None, y=None, data=None, width=500, height=500, p=None, **plot_kw):

    data['x'] = data[x].apply(lambda x: x if len(x) == 1 else tuple(x), axis=1)

    if p is None:
        p = figure(x_range=FactorRange(*data['x'].unique()), 
                   width=width, height=height,
                   min_border=100, x_axis_label=x[0], y_axis_label=y)

    p.circle(x=jitter('x', width=0.6, range=p.x_range), y=y, line_color='white',
             source=data, size=5, alpha=0.5, **plot_kw)

    return p

def pairplot(data=None, cols=None, hue=None, tooltips=None, width=400, height=400, output=None):
    if cols is None:
        cols = data.select_dtypes('number').columns

    comb_data = []
    for col1, col2 in combinations(cols, 2):
        tmp_data = data[[col1, col2]].assign(var1=col1, var2=col2)
        tmp_data.columns = ['value1', 'value2', 'var1', 'var2']

        if tooltips is not None:
            tmp_data = pd.concat([tmp_data, data[tooltips]], axis=1)
        if hue is not None:
            tmp_data[hue] = data[hue]
            
        comb_data.append(tmp_data)

    comb_data = pd.concat(comb_data)

    g = BokehFacetGrid(data=comb_data, hue=hue, row='var1', col='var2', tooltips=tooltips)
    g.map(scatter, x='value1', y='value2')

    if output is None:
        g.save('pairplot.html')


    
