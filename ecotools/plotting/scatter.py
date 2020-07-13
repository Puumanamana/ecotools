from itertools import combinations

import pandas as pd

from bokeh.plotting import figure
from bokeh.transform import jitter
from bokeh.models import FactorRange

from ecotools.parsing import parse_config
from ecotools.util import fit_ellipse, fit_hull, elt_or_nothing
from ecotools.plotting.grid import BokehFacetGrid

CFG = parse_config()['bokeh']

def scatter(x=None, y=None, data=None, width=500, height=500, p=None,
            ellipse=False, hull=False, s=5, alpha=1, line_color='black', hue_order=None,
            title=None, xlabel=None, ylabel=None,
            **plot_kw):

    if xlabel is None:
        xlabel = x[0]
    if ylabel is None:
        ylabel = y
    
    if p is None:
        p = figure(width=width, height=height, min_border=100, title=title,
                   x_axis_label=x[0], y_axis_label=y)    

    if (ellipse or hull) and 'legend_field' in plot_kw:
        hue = plot_kw['legend_field']
        for hue_i in data[hue].unique():
            data_subset = data[data[hue] == hue_i]
            if len(data_subset) < 3:
                continue

            if hull:
                hull_pts = fit_hull(data_subset[[x[0], y]].to_numpy())

                p.patch(*hull_pts.tolist(),
                        color=data_subset['color'].iloc[0],
                        fill_alpha=0.05)
            else:
                (mu, a, b, angle) = fit_ellipse(data_subset[[x[0], y]].to_numpy())
                p.ellipse(x=[mu[0]], y=[mu[1]], width=[2*a], height=[2*b], angle=angle,
                          color=data_subset['color'].iloc[0], fill_alpha=0.05)

    p.circle(x=x[0], y=y, line_color=line_color,
             source=data, size=s, alpha=alpha, **plot_kw)

    return p

def swarmplot(x=None, y=None, data=None, width=500, height=500, p=None, hue_order=None,
              jitter_width=0.6, s=5, alpha=0.5, line_color='black', **plot_kw):

    data['x'] = data[x].apply(lambda x: x if len(x) == 1 else tuple(x), axis=1)

    if p is None:
        p = figure(x_range=FactorRange(*data['x'].unique()), 
                   width=width, height=height,
                   min_border=100, x_axis_label=x[0], y_axis_label=y)

    p.circle(x=jitter('x', width=jitter_width, range=p.x_range), y=y, line_color=line_color,
             source=data, size=s, alpha=alpha, **plot_kw)

    return p


def lineplot(x=None, y=None, data=None, hue=None, width=800, height=800, p=None, hue_order=None,
              s=10, alpha=1, line_color='black', name=None, **plot_kw):

    if len(x) > 1:
        hue = x[-1]

    numeric_data = data.groupby(x)[y].agg(y='mean',
                                          inf=lambda x: max(x.min(), x.mean()-x.std()),
                                          sup=lambda x: min(x.max(), x.mean()+x.std()))

    hover_data = data.groupby(x).agg(elt_or_nothing)

    if 'color' in hover_data.columns:
        hover_data['color'] = hover_data.groupby(x[-1]).color.transform(lambda x: x.dropna().iloc[0])
    hover_data = hover_data.dropna(how='any', axis=1, thresh=numeric_data.count().min())

    line_data = pd.concat([numeric_data, hover_data], axis=1).reset_index()
    line_data['x_axis'] = line_data[x].apply(lambda x: x if len(x) == 1 else tuple(x), axis=1)

    if p is None:
        if hue is not None:
            x_ticks = FactorRange(*line_data.x_axis.tolist())
        else:
            x_ticks = line_data.x_axis.tolist()
        p = figure(x_range=x_ticks, x_axis_label=x[0], y_axis_label=y,
                   width=width, height=height, min_border=100)

    lines_y = line_data.groupby(x[-1])['y']
    lines_y = (lines_y.ffill() + lines_y.bfill()) / 2

    if hue is not None:
        lines_color = line_data.groupby(x[-1])['color'].first()
        lines_x = [[(x_i, hue_i) for x_i in line_data[x[0]].unique()] for hue_i in lines_color.index]
        lines_y = lines_y.groupby(line_data[x[-1]]).agg(list).tolist()
    else:
        lines_color = 'black'
        lines_x = [list(line_data[x[0]].unique())]
        lines_y = [lines_y.tolist()]

    p.segment('x_axis', 'inf', 'x_axis', 'y', line_color='black', line_alpha=0.8, source=line_data)
    p.segment('x_axis', 'y', 'x_axis', 'sup', line_color='black', line_alpha=0.8, source=line_data)

    p.multi_line(lines_x, lines_y, color=lines_color, line_width=2)

    p.circle(x='x_axis', y='y', line_color=line_color, size=s, alpha=alpha,
             source=line_data, **plot_kw)

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


    
