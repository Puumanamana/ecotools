from pathlib import Path
import os
import sys

import numpy as np

from ecotools.parser import parse_args, parse_config, get_strata
from ecotools.load_and_convert import convert, load_h5, save_h5
from ecotools.group_distances import get_normalized_distances, boxplot_with_factor_pairwise
from ecotools.correlations import pairwise_scatter
from ecotools.taxa_stackplot import stackplot
from ecotools.pathogens import get_pathogens
from ecotools.nmds import run_nmds, plot_nmds
from ecotools.subset import subset, subset_taxa

# functions defined in:
# https://github.com/hawaiidatascience/metaflowmics/blob/master/metaflowmics/scripts/bokeh_viz.py
sys.path.append(os.path.expanduser("~/projects/metaflowmics/metaflowmics/scripts"))
from bokeh_viz import preproc_for_heatmap, clustermap

CLADE_DIR = os.path.expanduser("~/databases/functional_clades")

def main():
    '''
    '''

    args = parse_args()

    if not args.command:
        sys.exit("No command selected")

    cfg = parse_config(args.cfg_file)

    strata = get_strata(cfg)

    for stratum in strata:

        print("{}: {}".format(args.command, stratum))

        suffix = ''
        if args.subset:
            suffix = '_subset'
        h5_file = Path("{}_{}{}.h5".format(cfg.get('misc', 'h5'), stratum, suffix))
        fig_dir = Path("{}_{}{}".format(cfg.get('misc', 'fig_dir'), stratum, suffix))
        table_dir = Path("{}_{}{}".format(cfg.get('misc', 'table_dir'), stratum, suffix))

        thresh = cfg.get('misc', 'otu_thresh')
        factors = cfg.get('metadata', 'factor_names').split(',')

        fig_dir.mkdir(exist_ok=True, parents=True)
        table_dir.mkdir(exist_ok=True, parents=True)

        if not h5_file.is_file():
            convert(cfg)

        (shared, taxonomy, alpha_div, metadata) = load_h5(h5_file)

        if shared.size == 0:
            sys.exit('No samples/OTUs found in abundance table')

        if args.command == 'distances':
            distances = get_normalized_distances(shared, metadata, args.factor)
            boxplot_with_factor_pairwise(distances, args.factor, fig_dir=fig_dir)

        elif args.command == 'nmds':
            components = run_nmds(shared, taxonomy, table_dir)
            plot_nmds(components, metadata,
                      color=args.color, row=args.row, col=args.col,
                      fig_dir=fig_dir)

        elif args.command == 'correlations':
            pairwise_scatter(shared, taxonomy, alpha_div, metadata, factor=args.factor, fig_dir=fig_dir)

        elif args.command == 'heatmap':
            metadata = metadata[factors]

            if args.factor:
                shared = shared.groupby(metadata[args.factor]).mean()
                metadata = metadata.groupby(args.factor).agg(lambda x: x[0] if len(x)==1 else float('nan'))
                metadata.dropna(how='all', axis=1, inplace=True)

                if metadata.size == 0:
                    metadata = None

            shared.index.name = 'Group'

            if args.pathogens:
                sp_file = Path(cfg.get('microbial', 'species'))
                pathogens = get_pathogens(sp_file)
                pathogenic_otus = np.intersect1d(pathogens.index, taxonomy.index)

                shared = shared[pathogenic_otus]
                taxonomy = taxonomy.merge(pathogens, left_index=True, right_index=True)

                title = "biclutered_pathogen_heatmap-by-{}".format(args.factor)

            (data, info) = preproc_for_heatmap(shared, taxonomy, meta=metadata)

            if not args.pathogens:
                title = "biclutered_heatmap_{}-by-genus_top-{}".format(args.factor, len(data.otu.unique()))

            clustermap(data, info, fig_dir=fig_dir, thresh=thresh, title=title)

        elif args.command == 'stackplot':

            if args.clade:
                clade_file = '{}/{}_blast/{}_genus.txt'.format(CLADE_DIR, args.clade, args.clade)
                taxa_subset = subset_taxa(taxonomy, clade_file)

                taxonomy = taxonomy[taxa_subset]
                shared = shared.loc[:, taxa_subset]

            stackplot(shared, taxonomy, rank=args.rank, factor=metadata[args.factor],
                      fig_dir=fig_dir, top=args.top, norm=args.norm, clade=args.clade)

        elif args.command == 'subset':

            datasets = subset(shared, taxonomy, metadata, alpha_div,
                              sample_file=args.sample_ids,
                              sample_cond=args.sample_cond,
                              taxa_file=args.taxa_ids,
                              taxa_cond=args.taxa_cond)

            output = '{}/{}_subset.h5'.format(h5_file.parent, h5_file.stem)
            
            save_h5(*datasets, output=output)
                
                
