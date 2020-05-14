import numpy as np
from bokeh.plotting import figure
from bokeh.transform import factor_cmap
from bokeh.models.ranges import FactorRange
from bokeh.models import Band, ColumnDataSource, HoverTool
from bokeh.layouts import gridplot

from ecotools.decorators import bokeh_legend_out, bokeh_save, bokeh_facets, bokeh_cmap
from ecotools.parsing import parse_config

CFG = parse_config()['bokeh']

def barplot(x=None, y=None, data=None, stacked=False, width=500, height=500,
            p=None, **plot_kw):

    grouped = data.groupby(x)
    x_values = list(grouped.groups.keys())

    if p is None:
        p = figure(x_range=FactorRange(*x_values), width=width, height=height,
                   x_axis_label=x[0], y_axis_label=y)

    data = {
        'x': x_values,
        'y': grouped[y].agg('mean'),
    }

    if len(x) > 1:
        data['color'] = grouped['color'].agg('first')
        data[x[1]] = grouped[x[1]].agg('first')

    p.vbar(x='x', top='y', width=0.9, source=data, **plot_kw)

    return p
    
def stackplot(x=None, y=None, data=None, width=500, height=500, norm=False,
              p=None, **plot_kw):

    plot_kw.pop('legend_field', None)
    plot_kw.pop('fill_color', None)    

    grouped = data.groupby(x)
    count_table = grouped[y].agg('mean').unstack()

    if norm:
        count_table = (count_table.T / count_table.sum(axis=1)).T
    
    x_values = count_table.index.tolist()
    colors = list(data['color'].unique())

    if p is None:
        p = figure(x_range=FactorRange(*x_values), width=width, height=height,
                   x_axis_label=x[0], y_axis_label=y)

    data = {col: count_table[col] for col in count_table.columns}
    data['x'] = x_values
    
    p.vbar_stack(count_table.columns.tolist(), x='x', width=0.9, color=colors,
                 legend_label=count_table.columns.tolist(),
                 source=ColumnDataSource(data), **plot_kw)

    return p
    


@bokeh_save
@bokeh_facets
@bokeh_legend_out
@bokeh_cmap
def stacked_barplot(metagenome, relative=False,**kwargs):

    metagenome = metagenome.copy()
    table = metagenome.abundance.data

    if relative:
        table = (table.T / table.sum(axis=1)).T

    tooltips = zip(table.columns, '@'+table.columns+'{0.00%}')

    factors = table.index.tolist()

    width = max(800, metagenome.n_samples()*50) + CFG['padding']
    height = 400 + CFG['padding']

    colors = [kwargs['cmap'][x] for x in table.columns]

    p = figure(plot_height=height, width=width, min_border=CFG['padding'], x_range=factors,
               tools=CFG['tools'], tooltips=list(tooltips))
    
    p.vbar_stack(table.columns, x=table.index.name, source=table.reset_index(),
                 width=.8, color=colors,
                 legend_label=table.columns.tolist())

    p.xaxis.major_label_orientation = 'vertical'
    
    return p


@bokeh_save
@bokeh_cmap
def permtest_barplot(metagenome, data, hue=None, strata=None, threshold=0.05, **kwarg):
    '''
    index = [strata, comparison]
    cols = [lvl1, lvl2, -log(p_adj), R2, n1, n2]
    '''

    hue_values = data[hue]

    data['neg_log_pval'] = -np.log10(data['pval_adj'].astype(float))
    data.set_index(strata + [hue], inplace=True)

    plot_data = dict(
        x=sorted(data.index.tolist()),
        pval=data['pval_adj'].tolist(),
        neg_log10_pval=data.neg_log_pval.tolist(),
        R2=data['R2'].tolist(),
        n1=data.n1.tolist(),
        n2=data.n2.tolist(),
    )
    band_data = dict(
        x=[-1, len(data) + len(hue_values) + 1],
        thresh=[-np.log10(threshold)]*2,
        upper=[data.neg_log_pval.max()*1.2]*2
    )
    
    plots = []

    index_cmap = factor_cmap('x', palette=list(kwarg['cmap'].values()),
                             factors=sorted(hue_values.unique()), start=1, end=2)

    tooltips = [('R2', '@R2'),
                ('pval_adj', '@pval'),
                ('#samples in {} 1'.format(hue), '@n1'),
                ('#samples in {} 2'.format(hue), '@n2')]

    titles = {
        'neg_log10_pval': 'Permanova test: -log10[pval_adj]',
        'R2': 'Permanova test: R2 score'
    }

    for metric, title in titles.items():
        p = figure(title=title, x_range=FactorRange(*plot_data['x']), tools=CFG['tools'][1:])
        p.vbar(x='x', top=metric, width=0.9, source=plot_data, name=metric,
               line_color="white", fill_color=index_cmap)

        if metric == 'neg_log10_pval':
            p.line(x='x', y='thresh', source=band_data,
                   line_color='grey', line_dash='dashed', line_width=1, line_alpha=0.5,
                   legend_label="{:.0%} significance threshold".format(threshold))

            band = Band(base='x', lower='thresh', upper='upper', level='underlay',
                        source=ColumnDataSource(band_data),
                        line_dash='dashed', fill_alpha=0.5, line_width=1, line_color='black')
            p.add_layout(band)
            p.legend.background_fill_alpha = 0
            p.legend.border_line_alpha = 0

        p.add_tools(HoverTool(names=[metric], tooltips=list(tooltips)))

        p.xaxis.major_label_orientation = "vertical"
        p.xaxis.axis_label_text_font_size = "10pt"
        plots.append(p)

        
    grid = gridplot(plots, ncols=1, plot_width=1300, plot_height=350)

    return grid
