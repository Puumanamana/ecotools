from itertools import combinations

from bokeh.plotting import figure

from ecotools.parsing import parse_config
# from ecotools.util import get_palette
from ecotools.decorators import bokeh_save, bokeh_facets, bokeh_legend_out, bokeh_cmap

CFG = parse_config()['bokeh']

@bokeh_save
@bokeh_facets
@bokeh_legend_out
@bokeh_cmap
def scatterplot(metagenome, values, hue=None, **kwargs):

    metadata =  metagenome.metadata.data[metagenome.metadata.qual_vars]
    data = values.merge(metadata, left_index=True, right_index=True).sort_values(by=hue)
    
    n_dims = values.shape[1]

    # hue_values = data[hue].unique()
    # palette = dict(zip(hue_values, get_palette(len(hue_values))))
    
    data['color'] = [kwargs['cmap'][x] for x in data[hue]]
    data.reset_index(inplace=True)

    tooltips = zip(data.drop('color', axis=1).columns, '@'+data.drop('color', axis=1).columns)

    plots = []
    for (k1, k2) in combinations(range(1, n_dims+1), 2):
        (label1, label2) = ('nmds_{}'.format(k1), 'nmds_{}'.format(k2))

        p = figure(title="(dim {} vs {})".format(k1, k2),
                   tooltips=list(tooltips), min_border=CFG['padding'])
        p.circle(x=label1, y=label2, color='color', line_color='gray',
                 source=data, size=10, alpha=0.5, legend_field=hue)
        plots.append(p)

    return plots
