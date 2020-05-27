from pathlib import Path

import pandas as pd

from bokeh.plotting import figure
from bokeh.models.ranges import FactorRange

from ecotools.plotting.facetgrid import BokehFacetGrid
from ecotools.parsing import parse_config
from ecotools.util import elt_or_nothing
from ecotools.decorators import timer

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
    p.xaxis.major_label_orientation = "vertical"

    return p


def stackplot(x=None, y=None, data=None, width=500, height=500, norm=False,
              hue_order=None, p=None,
              **plot_kw):

    for key in ['legend_field', 'fill_color']:
        plot_kw.pop('legend_field', None)

    # If x is provided, we aggregate
    agg_fn = dict([(col, 'mean') if col==y else (col, elt_or_nothing)
                   for col in data.columns if col not in x])

    x_sums = data.groupby(x[:-1])[y].sum()

    data = (data.groupby(x).agg(agg_fn)
            .dropna(subset=[y], how='any'))
    data.dropna(inplace=True, axis=1, thresh=data.count().max()*0.9)
    
    x_values = data.droplevel(x[-1]).index

    # Normalize to ratio to better visualize distribution
    if norm:
        data[y] = data.groupby(level=x[:-1])[y].transform(lambda x: x/x.sum())
        data['proportion'] = data[y].map(lambda x: '{:.1%}'.format(x))
    else:
        data['abundance'] = data[y].map(lambda x: '{:,}'.format(x))

    # Interactive features
    tooltips = list(zip(data.columns.drop([y, 'color']),
                        '@'+data.columns.drop([y, 'color'])))

    if x_sums.astype(float).max() > 1.1: # we have the abundance information
        tooltips = [('Sample size', '@group_size')] + tooltips

    # If it's the first plot to overlay
    if p is None:
        p = figure(x_range=FactorRange(*list(x_values.unique())),
                   width=width, height=height, min_border=100,
                   x_axis_label=x[0], y_axis_label=y, tooltips=tooltips)
    
    # Define all stacked bars
    table = data[y].unstack().fillna(0).sort_index(axis=1)

    data['bottom'] = table.shift(axis=1).cumsum(axis=1).fillna(0).stack()
    data['top'] = table.cumsum(axis=1).stack()

    data['group_size'] = (
        x_sums.map(lambda x: f'{x:,}').loc[x_values].to_numpy()
    )

    data = data.swaplevel(0, 1).sort_index()

    # plot each level one by one
    for level in data.index.get_level_values(x[-1]).unique():
        data_i = data.loc[level].assign(**{x[-1]: level})
        to_keep = data_i[['bottom', 'top']].dropna(how='any').index
        data_i = data_i.loc[to_keep]

        p.vbar(bottom='bottom',
               top='top',
               x=x[0], width=0.8, color='color',
               line_color='black', line_width=1.2,
               source=data_i.reset_index(),
               legend_label=level,
               name=level)

    p.xaxis.major_label_orientation = "vertical"

    return p
    

@timer
def taxa_stackplot(feature_table=None, feature_info=None, metagenome=None,
                   x='variable', row=None, col=None,
                   output='stackplot.html', norm=True, abd_thresh=0.01,
                   plot_kw={}, bar_kw={'norm': True}):
    
    '''Stacked barplot by sample groups
    
    Args:
        metagenome (MetagenomeDS): If the other dataframe information is skipped
        feature_table (pd.DataFrame): Count table (sample x OTU)
        taxonomy (pd.DataFrame): Taxomomy table (OTU x ranks)
        metadata (pd.DataFrame): Metadata table (sample x factors)
        norm (bool): Normalize sample group into ratios
        abd_thresh (float): Abundance threshold to group taxa into "others"
                            Must be in ]0, 1[
    Returns:
        None
    '''

    if metagenome:
        table = metagenome.get_column_format().reset_index()
        tax_cols = [table.columns[1]] + list(metagenome.taxonomy.columns)
        sample_var = table.columns[0]
    else:
        sample_var = 'variable'
        tax_cols = ['feature'] + list(feature_info.columns)
        table = (feature_table.T.merge(feature_info, left_index=True, right_index=True)
                 .rename_axis(index='feature').reset_index()
                 .melt(id_vars=tax_cols))

    hue = tax_cols[0]
    groups = [xi for xi in [x, row, col] if xi is not None]

    # Set threshold for assigning low abundance OTUs to others
    taxa_means = table.groupby(groups+[hue])['value'].agg('mean')
    sample_lims = taxa_means.sum(level=groups) * abd_thresh
    taxa_means = taxa_means.loc[table[groups+[hue]].apply(tuple, axis=1)]

    if len(groups) > 1:
        sample_lims = sample_lims.reindex(index=table[groups].apply(tuple, axis=1))
    else:
        sample_lims = sample_lims.reindex(index=table[groups[0]])

    in_others_cond = sample_lims.to_numpy() > taxa_means.to_numpy()

    table.loc[in_others_cond, tax_cols] = '_Others_'

    # Sum "others"
    agg_fn = dict((col, 'sum') if col == 'value' else (col, 'first')
                  for col in table.columns.drop([sample_var, hue]))
    table = table.groupby([sample_var, hue], as_index=False).agg(agg_fn)
    
    # Rank by total abundance
    hue_order = table.groupby(hue).value.sum().sort_values(ascending=False).index
    hue_order = hue_order.drop('_Others_').append(pd.Index(['_Others_']))

    g = BokehFacetGrid(data=table, row=row, col=col, hue=hue, hue_order=hue_order,
                       outdir=Path(output).parent, width=1000, **plot_kw)
    g.map(stackplot, x=x, y='value', **bar_kw)
    g.save(Path(output).name)

    
def stats_barplot(data, x=None, variables=['log10_p-adj', 'R2'], hue=None, threshold=0.05,
                  outdir='./', output='stats_barplot.html', plot_kw={}, bar_kw={}):

    data = data.melt(id_vars=[x for x in data.columns if x not in variables])

    g = BokehFacetGrid(data=data, hue=hue, row='variable', row_order=variables, outdir=outdir,
                       **plot_kw)
    g.map(barplot, x=x, y='value', tooltips=data.columns, **bar_kw)
    
    g.save(output)
    
