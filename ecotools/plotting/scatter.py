from itertools import combinations

import pandas as pd
from bokeh.plotting import figure
from bokeh.transform import jitter
from bokeh.models import FactorRange

from ecotools.parsing import parse_config
from ecotools.util import fit_ellipse
from ecotools.plotting.facetgrid import BokehFacetGrid

CFG = parse_config()['bokeh']

def scatter(x=None, y=None, hue=None, data=None, width=500, height=500, p=None,
            ellipse=False, s=5, alpha=.5, line_color='black', 
            **plot_kw):

    if hue is not None:
        plot_kw['legend_field'] = hue
        plot_kw['fill_color'] = 'color'

    if p is None:
        p = figure(width=width, height=height, min_border=100,
                   x_axis_label=x[0], y_axis_label=y)

    if ellipse:
        for hue_i in data[x[-1]].unique():
            data_subset = data[data[x[-1]] == hue_i]
            if len(data_subset) < 3:
                continue

            (mu, a, b, angle) = fit_ellipse(data_subset[[x[0], y]].to_numpy())
            p.ellipse(x=[mu[0]], y=[mu[1]], width=[2*a], height=[2*b], angle=angle,
                      fill_color=data_subset['color'].iloc[0],
                      line_color=line_color, alpha=alpha)

    p.circle(x=x[0], y=y, line_color=line_color,
             source=data, size=10, alpha=0.7, **plot_kw)

    return p

def swarmplot(x=None, y=None, data=None, width=500, height=500, p=None,
              jitter_width=0.6, s=5, alpha=0.5, line_color='black', **plot_kw):

    data['x'] = data[x].apply(lambda x: x if len(x) == 1 else tuple(x), axis=1)

    if p is None:
        p = figure(x_range=FactorRange(*data['x'].unique()), 
                   width=width, height=height,
                   min_border=100, x_axis_label=x[0], y_axis_label=y)

    p.circle(x=jitter('x', width=jitter_width, range=p.x_range), y=y, line_color=line_color,
             source=data, size=s, alpha=alpha, **plot_kw)

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


    
