from bokeh.plotting import figure

from ecotools.decorators import bokeh_legend_out, bokeh_save, bokeh_facets, bokeh_cmap
from ecotools.parsing import parse_config

CFG = parse_config()['bokeh']

@bokeh_save
@bokeh_facets
@bokeh_legend_out
@bokeh_cmap
def stacked_barplot(metagenome, **kwargs):

    metagenome = metagenome.copy()
    table = metagenome.abundance.data

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
