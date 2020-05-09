import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list

from bokeh.plotting import figure
from bokeh.models import ColorBar, LinearColorMapper, BasicTicker
from bokeh.transform import transform

from ecotools.parsing import parse_config
from ecotools.decorators import bokeh_save, bokeh_facets

CFG = parse_config()['bokeh']

def get_clusters(z_data, cluster_samples=True, cluster_features=True,
                 method='average', metric='braycurtis', optimal_ordering=True):

    sample_links = np.arange(z_data.shape[0])
    feature_links = np.arange(z_data.shape[1])

    if cluster_samples:
        try:
            sample_links = leaves_list(
                linkage(z_data, method=method, metric=metric, optimal_ordering=optimal_ordering)
            )
        except ValueError:
            print('Something went wrong with the hierarchical clustering on samples')

    if cluster_features:
        try:
            feature_links = leaves_list(
                linkage(z_data.T, method=method, metric=metric, optimal_ordering=optimal_ordering)
            )
        except ValueError:
            print('Something went wrong with the hierarchical clustering on features')

    return (z_data.index[sample_links], z_data.columns[feature_links])

@bokeh_save
@bokeh_facets
def clustermap(metagenome, table=None, groups=None, value_col='scores', standardize=True,
               cluster_samples=True, cluster_features=True, **kwargs):

    if table is None:
        table = metagenome.get_column_format()

    metagenome = metagenome.copy()

    if standardize:
        (samples, features) = table.index.names
        scores = table.groupby(features)[table.columns[0]].transform(lambda x: (x-x.mean()) / x.std())
    else:
        scores = table[value_col]

    (samples, features) = get_clusters(
        scores.unstack(), cluster_samples=cluster_samples, cluster_features=cluster_features
    )
    table['scores'] = scores

    kwargs.pop('output', None)

    p = heatmap(metagenome, table, 'scores', sample_order=samples, feature_order=features, **kwargs)
    return p


def heatmap(metagenome, data, value_col, sample_order=None, feature_order=None, **kwargs):

    (samples, features) = data.index.names

    if sample_order is None:
        sample_order = data.index.get_level_values(samples).unique()
    if feature_order is None:
        feature_order = data.index.get_level_values(features).unique()

    tooltips = dict(zip(data.columns, '@'+data.columns))

    p = figure(plot_height=max(500, 10*len(feature_order)),
               plot_width=max(500, len(sample_order)*10),
               x_range=sample_order.tolist(), y_range=feature_order.tolist(),
               x_axis_location="above",
               min_border=CFG['padding'],
               tools=CFG['tools'], tooltips=list(tooltips.items()))

    mapper = LinearColorMapper(palette='Magma256',
                               low=data[value_col].min(),
                               high=data[value_col].max())

    p.rect(x=samples, y=features, width=1, height=1, source=data.reset_index(),
           line_color=None, fill_color=transform(value_col, mapper))

    color_bar = ColorBar(color_mapper=mapper, formatter=p.xaxis.formatter,
                         ticker=BasicTicker(desired_num_ticks=5),
                         location=(0,0))
    
    p.add_layout(color_bar, 'right')

    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "8pt"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = 'vertical'

    return p
