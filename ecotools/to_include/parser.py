from pathlib import Path
from itertools import product
import os
import argparse
from configparser import ConfigParser, ExtendedInterpolation

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
    cfg.optionxform = str
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

def get_strata(cfg):
    strata_keys = set(cfg['strata'].keys()).difference(set(cfg['DEFAULT'].keys()))

    strata = {k: cfg.get('strata', k).split(',') for k in strata_keys}
    strata_prod = {'-'.join(comb): zip(strata.keys(), comb)
                   for comb in product(*strata.values())}

    return strata_prod

def parse_args():
    '''
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('--cfg-file', type=str, default=Path('../ikewai.ini'), action=ToPathAction, help='')
    parser.add_argument('--subset', action='store_true', help='Select last data subset')

    subparsers = parser.add_subparsers(help='sub-command help', dest='command')

    ss = subparsers.add_parser('subset', help='subset data based on indices or conditions')
    ss.add_argument('--taxa-ids', type=str, default=None, help='Path to column file with list of taxonomies. The column title must be the name of the rank to subset on')
    ss.add_argument('--sample-ids', type=str, default=None, help='Path to column file with list of samples. No header.')
    ss.add_argument('--taxa-cond', type=str, help='Taxa to keep based on a condition. Must be expressed as a list of {rank name}=={rank value}', nargs='*')
    ss.add_argument('--sample-cond', type=str, help='Sample to keep based on a condition. Must be expressed as a list of {metadata col}=={metadata value}', nargs='*')

    dist = subparsers.add_parser('distances', help='Plot distance distributions')
    dist.add_argument('--factor', type=str, default='Aquifer')

    nmds = subparsers.add_parser('nmds', help='Compute NMDS components')
    nmds.add_argument('--with-permanova', action='store_true', default=False)
    nmds.add_argument('--color', type=str, default='Aquifer', help='Color factor for scatter plot')
    nmds.add_argument('--row', type=str, help='Row factor for scatter plot')    
    nmds.add_argument('--col', type=str, help='Column factor for scatter plot')

    corr = subparsers.add_parser('correlations', help='Plot correlations between quantitative variables')
    corr.add_argument('--factor', type=str, default='Aquifer')
    
    hm = subparsers.add_parser('heatmap', help='Bi-clustered heatmap factor x OTUs')
    hm.add_argument('--factor', type=str, default='Aquifer')
    hm.add_argument('--pathogens', action='store_true', default=False)    

    sp = subparsers.add_parser('stackplot', help='Stacked barplot of OTU proportions per factor level')
    sp.add_argument('--factor', type=str, default='Aquifer')
    sp.add_argument('--rank', type=str, default='Phylum')
    sp.add_argument('--top', type=int, default=20)
    sp.add_argument('--norm', action='store_true', default=False)
    sp.add_argument('--clade', type=str, default=None, choices=['nirS', 'dsrA'])

    args = parser.parse_args()

    return args
