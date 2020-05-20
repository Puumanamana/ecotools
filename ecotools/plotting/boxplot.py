import pandas as pd

from bokeh.plotting import figure
from bokeh.models.ranges import FactorRange

from ecotools.parsing import parse_config
from ecotools.util import filter_metagenome, elt_or_nothing
from ecotools.plotting.facetgrid import BokehFacetGrid
from ecotools.plotting.scatter import swarmplot

CFG = parse_config()['bokeh']

def boxplot(x=None, y=None, data=None, width=500, height=500, p=None, **plot_kw):

    grouped = data.groupby(x)
    meta = grouped.agg(elt_or_nothing).dropna(axis=1)
    box_data = {col: meta[col] for col in meta.columns}

    quantiles = grouped[y].quantile(q=[0, 0.25, 0.5, 0.75, 1]).unstack()
    iqr = quantiles[0.75] - quantiles[0.25]

    box_data.update(dict(
        q1=quantiles[0.25], q2=quantiles[0.5], q3=quantiles[0.75], iqr=iqr,
        inf=grouped[y].agg(min), sup=grouped[y].agg(max),
        upper=quantiles[0.75] + 1.5*iqr, lower=quantiles[0.25] - 1.5*iqr,
        x=quantiles.index.tolist(),
    ))

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
    w_height = iqr.median() / 100
    p.rect('x', 'lower', 0.2, w_height, line_color="black", source=box_data)
    p.rect('x', 'upper', 0.2, w_height, line_color="black", source=box_data)

    p.xgrid.grid_line_color = None
    p.xaxis.major_label_orientation = "vertical"

    return p
    

def diversity_plot(metagenome, x=None, y=None, hue=None, col=None,
                   output='boxplot.html',
                   preproc_kw={}, plot_kw={}, box_kw={}, scatter_kw={}):

    mg = filter_metagenome(metagenome, inplace=False, **preproc_kw)
    mg.compute_alpha_diversity(y)
    data = pd.concat([mg.alpha_diversity, mg.metadata.factor_data()], axis=1)

    g = BokehFacetGrid(data=data, hue=hue, col=col, outdir=metagenome.figdir, **plot_kw)
    g.map(boxplot, x=x, y=y, **box_kw)
    g.map(swarmplot, x=x, y=y, tooltips=mg.metadata.qual_vars, **scatter_kw)
    g.save(output)




    