from pathlib import Path

from bokeh.plotting import figure
from bokeh.models.ranges import FactorRange

from ecotools.plotting.facetgrid import BokehFacetGrid
from ecotools.parsing import parse_config
from ecotools.util import elt_or_nothing

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
              p=None, **plot_kw):

    plot_kw.pop('legend_field', None)
    plot_kw.pop('fill_color', None)    

    # If x is provided, we aggregate
    agg_fn = dict([(col, 'mean') if col==y else (col, elt_or_nothing)
                   for col in data.columns if col not in x])
    data = data.groupby(x, as_index=False).agg(agg_fn).dropna(axis=1, how='any')

    # Normalize to ratio to better visualize distribution
    if norm:
        data[y] = data.groupby(x[:-1])[y].transform(lambda x: x/x.sum())

    x_values = list(data.groupby(x[:-1]).groups.keys())

    # Interactive features
    tips = data.drop(columns=x[:-1]+['color']).columns
    tooltips = list(zip(tips, '@'+tips))

    # If it's the first plot to overlay
    if p is None:
        p = figure(x_range=FactorRange(*x_values), width=width, height=height,
                   x_axis_label=x[0], y_axis_label=y, tooltips=tooltips)

    # Define all stacked bars
    grouped = data.groupby(x[0])[y]
    data['bottom'] = grouped.transform(lambda x: x.shift().cumsum()).fillna(0)
    data['top'] = grouped.cumsum()
    
    data_grouped = (data.groupby(x)
                    .agg('first')
                    .swaplevel(0, 1))

    stack_levels = list(data[x[-1]].unique())

    # plot each level one by one
    for i, level in enumerate(stack_levels):
        data_i = data_grouped.loc[level].assign(**{x[-1]: level})
        data_i['x'] = data_i.index

        p.vbar(bottom='bottom',
               top='top',
               x='x', width=0.8, color='color',
               source=data_i.reset_index(),
               legend_label=level,
               name=level)

    p.xaxis.major_label_orientation = "vertical"

    return p
    

def taxa_stackplot(feature_table=None, feature_info=None, metagenome=None,
                   x='variable', row=None, col=None,
                   output='stackplot.html', norm=True, abd_thresh=0.05,
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
        
    else:
        tax_cols = ['feature'] + list(feature_info.columns)
        table = (feature_table.T.merge(feature_info, left_index=True, right_index=True)
                 .rename_axis(index='feature').reset_index()
                 .melt(id_vars=tax_cols))

    grouping = [xi for xi in [x, row, col] if xi is not None]

    sample_lims = (table.groupby(grouping)['value'].agg(sum) * abd_thresh)

    if len(grouping) > 1:
        sample_lims = sample_lims.reindex(index=table[grouping].apply(tuple, axis=1))
    else:
        sample_lims = sample_lims.reindex(index=table[grouping[0]])

    table.loc[table.value < sample_lims.to_numpy(), tax_cols] = '_Others_'

    g = BokehFacetGrid(data=table, hue=tax_cols[0], row=row, col=col,
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
    
