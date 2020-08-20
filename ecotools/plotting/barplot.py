from pathlib import Path

import pandas as pd

from bokeh.plotting import figure
from bokeh.models.ranges import FactorRange

from ecotools.plotting.grid import BokehFacetGrid
from ecotools.parsing import parse_config
from ecotools.util import filter_groups

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
        data['color'] = grouped['color'].nth(0, dropna='all') 
        data[x[1]] = grouped[x[1]].nth(0, dropna='all') 

    p.vbar(x='x', top='y', width=0.9, source=data, **plot_kw)
    p.xaxis.major_label_orientation = "vertical"

    return p


def stackplot(x=None, y=None, data=None, width=1000, height=800, norm=True,
              hue_order=None, p=None,
              **plot_kw):

    for key in ['legend_field', 'fill_color']:
        plot_kw.pop('legend_field', None)

    x_sums = data.groupby(x[:-1])[y].sum()
    
    # Get uniq elt from metadata in each group, average value_var across x
    data = (data.groupby(x)
            .pipe(filter_groups, numeric=[y], fn='mean', approx=True)
            .dropna(subset=[y], how='any'))

    # Normalize to ratio to better visualize distribution
    if norm:
        data[y] = data.groupby(level=x[:-1])[y].transform(lambda x: x/x.sum())
        data['proportion'] = data[y].map(lambda x: '{:.1%}'.format(x))
    else:
        data['abundance'] = data[y].map(lambda x: '{:,}'.format(x))

    # Interactive features
    tooltips = list(zip(data.columns.drop([y, 'color', ' '], errors='ignore'),
                        '@'+data.columns.drop([y, 'color', ' '], errors='ignore')))

    if x_sums.astype(float).max() > 1.1: # we have the abundance information
        tooltips = [('Sample size', '@group_size')] + tooltips

    # Define all stacked bars
    table = data[y].unstack().fillna(0).sort_index(axis=1)

    data['bottom'] = table.shift(axis=1).cumsum(axis=1).fillna(0).stack()
    data['top'] = table.cumsum(axis=1).stack()

    # data = data.set_index(x).sort_index()
    # data = data.sort_values(by=x).set_index(x[-1])
    data['x_axis'] = data.reset_index(level=x[-1]).index

    data = data.assign(
        group_size=x_sums.map(lambda x: f'{x:,}').loc[data.x_axis].values,
        stack_name=data.index.get_level_values(x[-1])
    ).sort_index().reset_index(x[:-1])
    
    # data['x_axis'] = data.reset_index(level=x[-1]).index
    # data['group_size'] = x_sums.map(lambda x: f'{x:,}').loc[data.x_axis].values
    # data[x[-1]] = data.index
    # data = data.sort_index().reset_index(x[::-1])
    
    # If it's the first plot to overlay
    if p is None:
        p = figure(x_range=FactorRange(*list(data.x_axis.unique())),
                   width=width, height=height, min_border=100,
                   x_axis_label=x[0], y_axis_label=y, tooltips=tooltips)

    # plot each level one by one
    for level in data.stack_name.sort_values().unique():
        data_i = data.loc[[level]].dropna(how='any', subset=['bottom', 'top'])
        
        # data_i = data.xs(level, level=x[-1]).assign(**{x[-1]: level})
        # to_keep = data_i[['bottom', 'top']].dropna(how='any').index
        # data_i = data_i.loc[to_keep]
        # data_i['x'] = data_i.index

        p.vbar(bottom='bottom', top='top', x='x_axis',
               width=0.8, color='color',
               line_color='black', line_width=1.2,
               source=data_i,
               legend_label=level,
               name=level)
        
    p.xaxis.major_label_orientation = "vertical"

    return p
    

def taxa_stackplot(feature_table=None, feature_info=None, metagenome=None,
                   x='variable', hue=None, row=None, col=None,
                   output='stackplot.html', norm=True, abd_thresh=0.01,
                   plot_kw={}, bar_kw={}):
    
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

    x = [xi for xi in [x, hue] if xi is not None]

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
    groups = x + [xi for xi in [row, col] if xi is not None]

    # Set threshold for assigning low abundance OTUs to others
    taxa_means = table.groupby(groups+[hue])['value'].agg('mean')
    sample_lims = taxa_means.sum(level=groups) * abd_thresh
    taxa_means = taxa_means.loc[
        list(zip(*[table[x] for x in groups+[hue]]))
    ]
    
    if len(groups) > 1:
        sample_lims = sample_lims.reindex(index=table[groups].apply(tuple, axis=1))
    else:
        sample_lims = sample_lims.reindex(index=table[groups[0]])

    in_others_cond = sample_lims.to_numpy() > taxa_means.to_numpy()

    filler = 'Others (< {:.0%})'.format(abd_thresh)
    table.loc[in_others_cond, tax_cols] = filler
    
    agg_values = table.groupby([sample_var, hue]).value.sum()
    
    table = (table.groupby([sample_var, hue]).nth(0, dropna='all') # nth much faster than first()
             .assign(value=agg_values).reset_index())
    
    # Rank by total abundance
    hue_order = table.groupby(hue).value.sum().sort_values(ascending=False).index
    
    if filler in hue_order:
        hue_order = hue_order.drop(filler).append(pd.Index([filler]))

    g = BokehFacetGrid(data=table, row=row, col=col, hue=hue, hue_order=hue_order,
                       outdir=Path(output).parent, **plot_kw)
    g.map(stackplot, x=x, y='value', **bar_kw)
    g.save(Path(output).name)
    
def stats_barplot(data, x=None, variables=['log10_p-adj', 'R2'], hue=None, threshold=0.05,
                  outdir='./', output='stats_barplot.html', plot_kw={}, bar_kw={}):

    data = data.melt(id_vars=[x for x in data.columns if x not in variables])

    g = BokehFacetGrid(data=data, hue=hue, row='variable', row_order=variables, outdir=outdir,
                       **plot_kw)
    g.map(barplot, x=x, y='value', tooltips=data.columns, **bar_kw)
    
    g.save(output)
    
