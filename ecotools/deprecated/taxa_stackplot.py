from itertools import chain

from bokeh.io import output_file, save
from bokeh.plotting import figure
from bokeh.palettes import Set2, Set3, Category20, Accent, Dark2, Paired, viridis
from bokeh.models.tools import HoverTool

from ecotools.util import group_by_rank

TOOLS = ['hover', 'box_zoom']

def stackplot(shared, tax, factor, rank='Phylum', top=20, fig_dir=None, norm=False, clade=None):

    shared_group = shared.groupby(factor).agg('mean')
    tax_sub = tax.copy()
    
    if rank.lower().startswith('proteobacteria'):
        shared_group = shared_group.loc[:, tax.Phylum == 'Proteobacteria']
        tax_sub = tax.loc[shared_group.columns]
        rank = rank.split('_')[1]

    (shared_group, tax) = group_by_rank(shared_group, tax_sub, rank)

    if top is None:
        top = shared_group.shape[1]
        plot_name = "{}/barplot_by-{}_per-{}.html".format(fig_dir, rank, factor.name)
    else:
        top = min(top, shared_group.shape[1])
        plot_name = "{}/barplot_by-{}_per-{}_top-{}.html".format(fig_dir, rank, factor.name, top)

    criteria = shared_group.max().sort_values(ascending=False)
    shared_sub = shared_group.loc[:, criteria.index[:top]]

    col_order = shared_sub.sum().sort_values(ascending=False).index
    shared_sub = shared_sub.reindex(columns=col_order)
    
    if norm:
        shared_sub = (shared_sub.T/shared_sub.sum(axis=1)).T
        tooltips = zip(shared_sub.columns, '@'+shared_sub.columns+'{0.00%}')
        plot_name = plot_name.replace('.html', '_normed.html')

    if clade is not None:
        plot_name = plot_name.replace('.html', '_{}.html'.format(clade))

    palette = Set2[8] + Set3[12] + Category20[20]
    palette = list(chain(*[x[len(x)]
                           for x in [Set2, Set3, Category20, Accent, Dark2, Paired]]))

    if top > len(palette):
        palette = viridis(top)
    
    space = 400

    tooltips = zip(shared_sub.columns, '@'+shared_sub.columns)
    p = figure(plot_height=2*space+max(800, shared_sub.shape[0]*100),
               min_border=space,
               plot_width=1000+2*space,
               x_range=[0, 1.4*shared_sub.sum(axis=1).max()],
               y_range=shared_sub.index.tolist(),
               title="Mean {} {}abundance across {} (top {}/{})"
               .format(rank.lower(), 'relative '*norm, factor.name.lower(), top, shared_group.shape[1]),
               tooltips=list(tooltips))

    p.hbar_stack(shared_sub.columns, y=factor.name, source=shared_sub.reset_index(),
                 height=.5, color=palette[:top],
                 legend_label=shared_sub.columns.tolist())

    p.select(dict(type=HoverTool)).point_policy = "follow_mouse"
    
    output_file(plot_name)
    save(p)

def elem_or_nothing(x):
    if len(x) == 1:
        return x[0]    
    return float('nan')
