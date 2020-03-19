import os
import sys
import argparse
from itertools import chain

from bokeh.io import output_file, save
from bokeh.plotting import figure
from bokeh.palettes import Set2, Set3, Category20, Accent, Dark2, Paired, viridis
from bokeh.models.tools import HoverTool

from load_and_convert import load_h5, parse_config

sys.path.append(os.path.expanduser("~/projects/metaflowmics/metaflowmics/scripts"))
from bokeh_viz import preproc_for_heatmap, clustermap

TOOLS = ['hover', 'box_zoom']
OUTDIR = '../outputs'

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--factor', type=str, default='Aquifer')
    parser.add_argument('--norm', action='store_true', default=False)
    args = parser.parse_args()

    return args

def stackplot(shared, tax, factor, rank='Phylum', top=20, fig_dir=None, norm=False):

    shared_sub = shared.groupby(factor).agg('mean')
    tax_sub = tax.copy()
    
    if rank == 'Class':
        shared_sub = shared_sub.loc[:, tax.Phylum == 'Proteobacteria']
        tax_sub = tax.loc[shared_sub.columns]

    shared_sub = shared_sub.T.groupby(tax_sub[rank]).agg('sum').T

    if top is None:
        top = shared_sub.shape[1]
        plot_name = "{}/barplot_by-{}_per-{}.html".format(fig_dir, rank, factor.name)
    else:
        top = min(top, shared_sub.shape[1])
        plot_name = "{}/barplot_by-{}_per-{}_top-{}.html".format(fig_dir, rank, factor.name, top)

    criteria = shared_sub.max().sort_values(ascending=False)
    shared_sub = shared_sub.loc[:, criteria.index[:top]]

    col_order = shared_sub.sum().sort_values(ascending=False).index
    shared_sub = shared_sub.reindex(columns=col_order)
    
    if norm:
        shared_sub = (shared_sub.T/shared_sub.sum(axis=1)).T
        tooltips = zip(shared_sub.columns, '@'+shared_sub.columns+'{0.00%}')
        plot_name = plot_name.replace('.html', '_normed.html')

    palette = Set2[8] + Set3[12] + Category20[20]
    palette = list(chain(*[x[len(x)] for x in [Set2, Set3, Category20, Accent, Dark2, Paired]]))

    if top > len(palette):
        palette = viridis(top)
    
    space = 400

    tooltips = zip(shared_sub.columns, '@'+shared_sub.columns)
    p = figure(plot_height=2*space+max(800, shared_sub.shape[0]*100),
               min_border=space,
               plot_width=1000+2*space,
               x_range=[0, 1.4*shared_sub.sum(axis=1).max()],
               y_range=shared_sub.index.tolist(),
               title="Mean {} {}abundance across {}"
               .format(rank.lower(), 'relative '*norm, factor.name.lower()),
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

if __name__ == '__main__':
    args = parse_args()
    cfg = parse_config()

    thresh = cfg.get('misc', 'otu_thresh')
    fig_dir = cfg.get('misc', 'fig_dir')
    
    (shared, taxonomy, _, metadata) = load_h5(cfg.get('misc', 'h5'))
    print('Data loaded')

    stackplot(shared, taxonomy, rank='Phylum', factor=metadata[args.factor],
              fig_dir=fig_dir, top=30, norm=args.norm)
    stackplot(shared, taxonomy, rank='Class', factor=metadata[args.factor], fig_dir=fig_dir, norm=args.norm)

    cols = cfg.get('metadata', 'factor_names').split(',')
    shared = shared.groupby(metadata[args.factor]).mean()
    shared.index.name = 'Group'

    if args.factor:
        metadata = metadata[cols].groupby(args.factor).agg(elem_or_nothing).dropna(how='all', axis=1)
        if metadata.size == 0:
            metadata = None
    else:
        metadata = metadata[cols]

    (data, info) = preproc_for_heatmap(shared, taxonomy, meta=metadata)
    title = "biclutered_heatmap_{}-by-genus-{}_top-{}".format(args.factor, thresh, len(data.otu.unique()))
    clustermap(data, info, fig_dir=fig_dir, thresh=thresh, title=title)
