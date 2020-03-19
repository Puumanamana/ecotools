from pathlib import Path
import os
import sys
import argparse
from configparser import ConfigParser, ExtendedInterpolation

from src.load_and_convert import convert, load_h5
from src.group_distances import get_normalized_distances, boxplot_with_factor_pairwise
from src.correlations import pairwise_scatter
from src.factors_vs_16S import stackplot

sys.path.append(os.path.expanduser("~/projects/metaflowmics/metaflowmics/scripts"))
from bokeh_viz import preproc_for_heatmap, clustermap

class ToPathAction(argparse.Action):
    '''
    argparse action to convert string to Path objects
    '''
    def __init__(self, option_strings, dest, required=False, **kwargs):
        argparse.Action.__init__(self, option_strings=option_strings, dest=dest, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        if isinstance(values, (list, tuple)):
            values_ok = [Path(val) for val in values]
        else:
            values_ok = Path(values)

        setattr(namespace, self.dest, values_ok)

def parse_config(filename):
    cfg = ConfigParser({'home': os.environ['HOME']},
                       interpolation=ExtendedInterpolation())
    cfg.read(filename)

    # Save parsed file for R
    with open('{}/R_{}'.format(filename.parent, filename.name), 'w') as configfile:
        default_fields = dict(cfg.items('DEFAULT')).keys()
        for section in cfg.sections():
            for (name, value) in cfg.items(section):
                if name not in default_fields:
                    cfg.set(section, name, value)
        cfg.write(configfile)

    return cfg    

def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--cfg-file', type=str, default=Path('../ikewai.ini'), action=ToPathAction,
                        help='')

    subparsers = parser.add_subparsers(help='sub-command help', dest='command')

    # subparsers.add_parser('convert',
    #                       help='Convert to h5 format')

    dist = subparsers.add_parser('distances',
                                 help='Plot distance distributions')
    dist.add_argument('--factor', type=str, default='Aquifer')

    corr = subparsers.add_parser('correlations',
                                 help='Plot correlations between quantitative variables')
    corr.add_argument('--factor', type=str, default='Aquifer')
    
    hm = subparsers.add_parser('heatmap',
                                 help='Bi-clustered heatmap factor x OTUs')
    hm.add_argument('--factor', type=str, default='Aquifer')

    sp = subparsers.add_parser('stackplot',
                               help='Stacked barplot of OTU proportions per factor level')
    sp.add_argument('--factor', type=str, default='Aquifer')
    sp.add_argument('--rank', type=str, default='Phylum')
    sp.add_argument('--top', type=int, default=20)
    sp.add_argument('--norm', action='store_true', default=False)

    args = parser.parse_args()

    return args

def main():
    '''
    '''

    args = parse_args()
    cfg = parse_config(args.cfg_file)

    h5_file = Path(cfg.get('misc', 'h5'))
    fig_dir = cfg.get('misc', 'fig_dir')
    thresh = cfg.get('misc', 'otu_thresh')
    factors = cfg.get('metadata', 'factor_names').split(',')
    

    if not h5_file.is_file():
        convert(cfg)
        
    (shared, taxonomy, alpha_div, metadata) = load_h5(h5_file)

    if args.command == 'distances':
        distances = get_normalized_distances(h5_file, args.factor)
        boxplot_with_factor_pairwise(distances, args.factor, fig_dir=fig_dir)

    elif args.command == 'correlations':
        pairwise_scatter(shared, taxonomy, alpha_div, metadata, factor=args.factor, fig_dir=fig_dir)

    elif args.command == 'heatmap':
        # following functions defined in
        # https://github.com/hawaiidatascience/metaflowmics/blob/master/metaflowmics/scripts/bokeh_viz.py
        if args.factor:
            shared = shared.groupby(metadata[args.factor]).mean()
            metadata = metadata.groupby(args.factor).agg(lambda x: x[0] if len(x)==1 else float('nan'))
            metadata.dropna(how='all', axis=1, inplace=True)

            if metadata.size == 0:
                metadata = None

        shared.index.name = 'Group'

        (data, info) = preproc_for_heatmap(shared, taxonomy, meta=metadata[factors])
        title = "biclutered_heatmap_{}-by-genus_top-{}".format(args.factor, len(data.otu.unique()))
        clustermap(data, info, fig_dir=fig_dir, thresh=thresh, title=title)

    elif args.command == 'stackplot':
        stackplot(shared, taxonomy, rank=args.rank, factor=metadata[args.factor],
                  fig_dir=fig_dir, top=args.top, norm=args.norm)

    elif args.command == 'metadata':
        pass
    

if __name__ == '__main__':
    main()
