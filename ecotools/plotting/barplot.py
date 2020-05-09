import numpy as np
from bokeh.plotting import figure
from bokeh.transform import factor_cmap
from bokeh.models.ranges import FactorRange
from bokeh.models import Band, ColumnDataSource
from bokeh.layouts import gridplot

from ecotools.decorators import bokeh_legend_out, bokeh_save, bokeh_facets, bokeh_cmap
from ecotools.parsing import parse_config

CFG = parse_config()['bokeh']

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
    
    data['neg_log_pval'] = -np.log10(data['p-value'].astype(float))
    data.set_index(strata + [hue], inplace=True)
    
    plot_data = dict(
        x=sorted(data.index.tolist()),
        pval=data['p-value'].tolist(),
        neg_log10_pval=data.neg_log_pval.tolist(),
        # R2=data['R2'].tolist(),
        upper=[data.neg_log_pval.max()*1.2]*len(data),
        lower=[-np.log10(threshold)]*len(data),
        n1=data.n1.tolist(),
        n2=data.n2.tolist(),
    )
    plots = []

    index_cmap = factor_cmap('x', palette=list(kwarg['cmap'].values()),
                             factors=sorted(hue_values.unique()), start=1, end=2)

    tooltips = [# ('R2', '@R2'),
                ('p-value', '@pval'),
                ('#samples in {} 1'.format(hue), '@n1'),
                ('#samples in {} 2'.format(hue), '@n2')]

    titles = {
        'neg_log10_pval': 'Permanova test: -log10[p-value]',
        # 'R2': 'Permanova test: R2 score'
    }
    
    for metric, title in titles.items():
        p = figure(x_range=FactorRange(*plot_data['x']), tooltips=tooltips)
        p.vbar(x='x', top=metric, width=0.9, source=plot_data,
               line_color="white", fill_color=index_cmap)

        if metric == 'neg_log10_pval':
            p.line(np.arange(len(data) + len(hue_values)),
                   [-np.log10(threshold)]*(len(data) + len(hue_values)),
                   line_color='grey', line_dash='dashed', line_width=1, line_alpha=0.5,
                   legend_label="{:.0%} significance threshold".format(threshold))

            band = Band(base='x', lower='lower', upper='upper', level='underlay',
                        source=ColumnDataSource(plot_data),
                        line_dash='dashed', fill_alpha=0.5, line_width=1, line_color='black')
            p.add_layout(band)
            p.legend.background_fill_alpha = 0.0

        p.xaxis.major_label_orientation = "vertical"
        p.xaxis.axis_label_text_font_size = "10pt"
        plots.append(p)

    grid = gridplot(plots, ncols=1, plot_width=1500, plot_height=500)

    return grid
