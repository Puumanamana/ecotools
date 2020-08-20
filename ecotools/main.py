import warnings
import sys
from pathlib import Path
import pandas as pd
import numpy as np

import ecotools
from ecotools.parsing import parse_args

from ecotools.plotting.boxplot import diversity_plot, boxplot
from ecotools.plotting.barplot import taxa_stackplot
from ecotools.plotting.heatmap import clustermap
from ecotools.plotting.grid import BokehFacetGrid
from ecotools.plotting.scatter import scatter, swarmplot

from ecotools.machine_learning.lda import lda_model, lda_boxplot
from ecotools.beta_diversity.ordination import ordinate

def load_metagenome(input_dir, metadata_path, output_dir, otu_thresh=100, metadata_parser={}):
    '''
    Args:
        input_dir (str): MetaFlowmics output directory.
        metadata_path (str): Path to metadata file. Must be csv formatted with header and rownames as unique sample id.
        otu_thresh (int): OTU clustering threshold. 100 is no clustering (ASVs)
        metadata_parser (dict): extra arguments for parsing metadata file. See ecotools.metagenome_from_files() arguments for more details

    Returns:
        MetagenomicsDS: metagenome object from ecotools

    '''     
    pipeline_files = ecotools.util.find_pipeline_files(input_dir)
    pipeline_files.pop('fasta_path', None)
    pipeline_files.pop('tree_path', None)    

    mg = ecotools.metagenome_from_files(
        **pipeline_files, metadata_path=metadata_path,
        outdir=output_dir, parser_args=metadata_parser
    )
    mg.to_h5()

    return mg


def preprocess_metagenome(mg, min_prevalence=0, subsample=False, relabund=False,
                          norm_func=None, otu_subset=None, inplace=True):
    '''
    Args:
        mg (MetagenomicsDS): metagenome
        min_prevalence (int): Minimum number of sample a given OTU must occur in
        subsample (bool): Subsample the counts. If True, the level is automatically determined
        relabund (bool): Convert counts to relative abundance
        norm_func (function): Convert raw counts using the nom_func function
        otu_subset (str): Path to taxonomy file to subset OTUs or [(rank, names) list]

    Returns:
        MetagenomicsDS
    '''
    if not inplace:
        mg = mg.copy()

    mg.subset_otus((mg.abundance.data > 0).sum() > min_prevalence)

    if subsample:
        mg.subsample()

    if relabund:
        mg.to_relative_abundance()

    if callable(norm_func):
        mg.abundance.data = norm_func(mg.abundance.data)

    if otu_subset and otu_subset is not None:
        if isinstance(otu_subset, str) or isinstance(otu_subset, Path):
            if Path(otu_subset).is_file():
                mg.subset_otus(otus=Path(otu_subset))
            else:
                print(f'{otu_subset} does not exist. Continuing without clade subset')
        elif isinstance(otu_subset, list):
            mg.subset_otus(taxa=otu_subset)
        else:
            print('Unknown subsetting object.')
            import ipdb;ipdb.set_trace()
            exit(1)
    
    print(mg)

    if not inplace:
        return mg


def distr_cmd(mg, factors, diversity=[], otu_list=[]):
    (hue, x, row, col) = factors

    if diversity:
        for div in diversity:
            print(f'Diversity boxplot: {div}')
            diversity_plot(mg, x=x, y=div, hue=hue, points=False,
                           col=col, row=row,
                           plot_kw={'height': 400},
                           output=f'boxplot_{div}.html')

    if otu_list:
        if len(otu_list) == 1 and Path(otu_list[0]).is_file():
            otu_list = open(otu_list[0]).read().splitlines()[1:]

        otu_list = list(set(otu_list).intersection(set(mg.columns)))[:50]

        if not otu_list:
            return
        
        print(f'Feature boxplot ({len(otu_list)} otus (max=50))')
        
        data = (mg.get_column_format(tax=False)
                .loc[(slice(None), otu_list), :]
                .reset_index())

        repl = {}
        for otu in otu_list:
            lineage = mg.taxonomy.data.loc[otu].tolist()[:6]
            short = [x[0].lower() for x in mg.taxonomy.ranks][:len(lineage)]
            new_name = ', '.join(['{}: {}'.format(*it) for it in zip(short, lineage)])
            repl[otu] = new_name

        data = data.replace(repl)

        g = BokehFacetGrid(data=data, hue=hue, col='OTU', col_wrap=5, outdir=mg.figdir,
                           width='auto', height=400)
        g.map(boxplot, x=x, y='value', tooltips=['group_size', 'not_null'])
        g.map(swarmplot, x=x, y='value', tooltips=data.columns)#.drop(columns=['OTU', 'color']))
        g.save('specific_otus_distribution.html')


def tax_cmd(mg, factors=None, ranks=None, plot_bars=True, plot_heatmap=True):
    (hue, x, row, col) = factors
    
    for rank in ranks:
        print(f'Plotting taxonomic composition for: {rank}')
        if plot_bars:
            mg_ = mg.copy()
            mg_.group_taxa(rank)

            taxa_stackplot(metagenome=mg_, x=x, hue=hue, col=col, row=row, norm=True,
                           output=f'{mg.figdir}/barplot-{rank}.html',
                           plot_kw={'width': 1400, 'height': 800})

        if plot_heatmap:
            mg_ = mg.copy()
            # mg_.abundance.data = np.sqrt(mg_.abundance.data)
            data = mg_.get_column_format()
            g = BokehFacetGrid(data=data, row=x, col=col, outdir=mg_.figdir)
            g.map(clustermap, y=hue, x=rank, z='value', cluster_samples=False)
            g.save(f'clustermap-{rank}.html')

def lda_cmd(mg, factors, k=10):
    lda_results = lda_model(mg, k=k)

    taxa_stackplot(output=f'{mg.figdir}/topic-feature_stackplot.html',
                   feature_table=lda_results['features'], abd_thresh=0.01,
                   feature_info=mg.taxonomy.data, plot_kw={'width': 1200, 'height': 800})

    (x, row, col, _) = factors
    lda_boxplot(data=lda_results,
                metadata=mg.metadata.factor_data(),
                taxonomy=mg.taxonomy.data.assign(OTU=mg.taxonomy.index),
                x=x, row=row, rank='OTU', top=10,
                output=f'{mg.figdir}/topic-sample_boxplot.html')

def ordination_cmd(mg, factors=[], strata=[], method='pcoa', distance='bray'):

    (hue, col, row, extra) = factors[:1] + strata + [None] * (4-len(strata)-1)

    if extra is not None or len(factors) > 1:
        extra = ','.join([str(extra)] + [str(x) for x in factors[1:]])
        warnings.warn(f'Could not render {extra} too many levels to render. Ignoring.', UserWarning)

    result = ordinate(mg, strata=strata, subsample=True)
    if not all(x is None for x in strata):
        components = pd.concat(val['sample'] for val in result.values())
    else:
        components = result['sample']

    compo_names = components.columns[:2]
    other_meta = np.setdiff1d(mg.metadata.qual_vars, components.columns)
    components = pd.concat([components, mg.metadata.factor_data(other_meta)], axis=1).dropna(subset=compo_names, how='any')

    hull = (len(components[hue].unique()) < 5)
    
    g = BokehFacetGrid(data=components.reset_index(), hue=hue, col=col, row=row, outdir=mg.figdir, scale=1.5)
    g.map(scatter, x=components.columns[0], y=components.columns[1], hull=hull, s=10,
          tooltips=mg.metadata.qual_vars)
    g.save(f'{method}-{distance}.html')

    g = BokehFacetGrid(
        data=components.melt(id_vars=components.columns.drop(compo_names),
                             value_name='score', var_name='component'),
        col=col, row=row, outdir=mg.figdir, hue=hue)
    g.map(boxplot, x='component', y='score')
    g.map(swarmplot, x='component', y='score', tooltips=components.columns[2:])    
    g.save(f'{method}_{distance}_boxplot.html')
    
    
def main():

    args = parse_args()

    mg = load_metagenome(args.input_dir, args.metadata, args.output,
                         otu_thresh=args.otu_thresh,
                         metadata_parser=dict(qual_vars=args.qual))

    preprocess_metagenome(mg, min_prevalence=args.min_prevalence,
                          subsample=args.subsample,
                          relabund=args.relabund,
                          otu_subset=args.otu_subset)
    
    if args.cmd.startswith('distr'):
        distr_cmd(mg, args.conditions, args.diversity, args.otu_list)
    elif args.cmd.startswith('taxonomy'):
        tax_cmd(mg, args.conditions, ranks=args.ranks,
                plot_bars=args.bar, plot_heatmap=args.heatmap)
    elif args.cmd.startswith('lda'):
        lda_cmd(mg, args.conditions, k=args.n_topics)
    elif args.cmd.startswith('ordination'):
        mg.run = dict(metric=args.distance, ordination=args.method, clade=args.otu_subset,
                      factors=args.conditions, strata=args.strata)
        ordination_cmd(mg, factors=args.conditions, strata=args.strata,
                       method=args.method, distance=args.distance)
    else:
        sys.exit(f'Unknown command: {args.cmd}')
    

if __name__ == '__main__':
    main()
