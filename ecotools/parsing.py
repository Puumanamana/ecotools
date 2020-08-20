from pathlib import Path
import yaml
import argparse

from ecotools.core.taxonomy_table import TaxonomyTable

SRC_DIR = Path(__file__).resolve().parent
TOOLS = ['taxonomy', 'distr', 'lda', 'ordination']

def parse_config():
    with open('{}/settings.yaml'.format(SRC_DIR)) as handle:
        cfg = yaml.load(handle, Loader=yaml.FullLoader)

    return cfg


def parse_args():

    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers(title='commands', dest='cmd')
    subparsers.required = True

    handles = {}

    for tool in TOOLS:
        p = subparsers.add_parser(tool)
        p.add_argument('-i', '--input-dir', help='Path to MetaFlowmics output folder')
        p.add_argument('-o', '--output', help='Path to output folder', default='ecotools_outputs')
        p.add_argument('-m', '--metadata', help='Path to metadata file (csv formatted)')
        p.add_argument('--qual', type=str, nargs='+', help='Int covariates to keep as factors')
        p.add_argument('--conditions', type=str, nargs='*', required=True)
        p.add_argument('--otu-subset', type=str, default ='', help='Path to clade file')
        p.add_argument('--otu-thresh', type=str, default=100)
        p.add_argument('--min-prevalence', type=int, default=0)
        p.add_argument('--subsample', action='store_true')
        p.add_argument('--relabund', action='store_true')
        handles[tool] = p

    handles['distr'].add_argument('--diversity', type=str, nargs='+', help='Diversity metrics to plot')
    handles['distr'].add_argument('--otu-list', type=str, nargs='*', help='Specific OTU distributions to plot')
    handles['taxonomy'].add_argument('--ranks', type=str, nargs='+', choices=TaxonomyTable.ranks, required=True)
    handles['taxonomy'].add_argument('--bar', action='store_true')
    handles['taxonomy'].add_argument('--heatmap', action='store_true')    
    handles['lda'].add_argument('--n-topics', type=int, default=10, help='Number of topics for LDA analysis')
    handles['ordination'].add_argument('--method', type=str, default='pcoa')
    handles['ordination'].add_argument('--distance', type=str, default='bray')
    handles['ordination'].add_argument('--strata', type=str, nargs='*', default=[])    
    
    args = parser.parse_args()

    args.conditions = args.conditions + [None] * (4-len(args.conditions))

    return args
