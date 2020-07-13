import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list

from bokeh.plotting import figure
from bokeh.models import ColorBar, LinearColorMapper, BasicTicker
from bokeh.transform import transform

from ecotools.parsing import parse_config
from ecotools.util import filter_metagenome
from ecotools.plotting.grid import BokehFacetGrid

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



def clustermap(x=None, y=None, z=None, data=None, standardize=True,
               cluster_samples=True, cluster_features=True, **plot_kw):

    data_uniq = data.groupby([x[0], y]).agg('first')
    data_uniq[z] = data.groupby([x[0], y])[z].agg('mean')
    data_uniq.reset_index(inplace=True)
    
    if standardize:
        scores = data_uniq.groupby(x[0])[z].transform(lambda x: (x-x.mean()) / x.std())
    else:
        scores = data_uniq[z]

    data_uniq['scores'] = scores
    data_uniq = data_uniq[data_uniq.scores.notnull()]

    (samples, features) = get_clusters(
        data_uniq[x + [y, z]].pivot(y, x[0], values=z),
        cluster_samples=cluster_samples, cluster_features=cluster_features
    )

    p = heatmap(data=data_uniq, x=x, y=y, z='scores', x_order=features, y_order=samples, **plot_kw)
    return p


def heatmap(x=None, y=None, z=None, data=None, x_order=None, y_order=None, p=None,
            width=1200, height=500, **plot_kw):

    if x_order is None:
        x_order = data[x[0]].unique()
    if y_order is None:
        y_order = data[y].unique()

    tooltips = zip(data.columns, '@'+data.columns)
        
    if p is None:
        p = figure(height=max(height, 10*len(y_order)),
                   plot_width=max(width, 10*len(x_order)),
                   tooltips=list(tooltips), tools=CFG['tools'],
                   x_range=x_order.tolist(), y_range=y_order.tolist(),
                   x_axis_location="above",
                   min_border=CFG['padding'])

    mapper = LinearColorMapper(palette='Magma256',
                               low=data[z].min(),
                               high=data[z].max())

    p.rect(x=x[0], y=y, width=1, height=1, source=data,
           line_color=None, fill_color=transform(z, mapper))

    color_bar = ColorBar(color_mapper=mapper,
                         ticker=BasicTicker(desired_num_ticks=5),
                         location=(0,0))
    
    p.add_layout(color_bar, 'right')

    p.axis.axis_line_color = None
    p.axis.major_tick_line_color = None
    p.axis.major_label_text_font_size = "8pt"
    p.axis.major_label_standoff = 0
    p.xaxis.major_label_orientation = 'vertical'

    return p


def metagenome_heatmap(mg, y=None, col=None, row=None, output='heatmap.html',
                       preproc_kw={}, **kwargs):

    mg = filter_metagenome(mg, inplace=False, **preproc_kw)        
    mg.taxonomy.clean_labels(trim=True)

    if y is None:
        y = 'group'

    x = 'OTU'
    if 'rank' in preproc_kw:
        x = preproc_kw['rank']

    data = mg.get_column_format().reset_index()

    g = BokehFacetGrid(data=data, outdir=mg.figdir, col=col, row=row)
    g.map(clustermap, x=x, y=y, z='value',
          standardize=True, cluster_rows=True, cluster_cols=True, **kwargs)
    g.save(output)

    
    
