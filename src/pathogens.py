from pathlib import Path

import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list

from load_and_convert import load_h5, parse_config
from factors_vs_16S import clustermap, parse_args

def get_pathogens(sp_file, db_file):
    species = pd.read_csv(sp_file, index_col=0).Species.dropna()
    db = pd.read_csv(db_file, header=None, sep='\t')[0].drop_duplicates()

    patho_otus = set()
    for i, sp in enumerate(species.to_numpy()):
        res = np.mean([db.str.contains(hit, case=False).sum() > 0 for hit in sp.split('/')])

        if res > 0.5:
            patho_otus.add(i)
            print(len(patho_otus), end='\r')

    print("{} / {} potential pathogenic OTUs".format(len(patho_otus), len(species)))
    return species.iloc[list(patho_otus)].sort_index()

def get_clustered_zscores(df, factor=None, metadata=None):

    if factor:
        df = df.groupby(metadata[factor]).agg('mean')
    else:
        factor = 'sampleID'

    df = (df - df.mean()) / df.std()
    df.dropna(inplace=True, axis=1)
    df.index.name = 'variable'

    # hierarchical clustering on both rows and columns
    try:
        row_links = leaves_list(linkage(df, method='average', metric='braycurtis', optimal_ordering=True))
    except ValueError:
        print('Something went wrong with the hierarchical clustering on rows')
        row_links = np.arange(df.shape[0])
    try:
        col_links = leaves_list(linkage(df.T, method='average', metric='braycurtis', optimal_ordering=True))
    except ValueError:
        print('Something went wrong with the hierarchical clustering on cols')
        col_links = np.arange(df.shape[1])

    df = (df.iloc[row_links, col_links].T
            .reset_index()
            .melt(id_vars=['index'])
            .rename(columns={'value': 'z_score', 'variable': factor, 'index': 'otu'})
    )

    return df

def get_shared_info(data, factor=None, metadata=None):
    abd = data.reset_index().melt(id_vars='index')

    if factor:
        abd = (abd.merge(metadata[factor], left_on='index', right_index=True)
               .groupby([factor, 'variable']).sum()
               .reset_index())
        abd.columns = [factor, 'otu', 'abundance']
        metadata = pd.DataFrame(metadata[factor].value_counts().rename('n_samples'))

    else:
        abd.columns = ['sampleID', 'otu', 'abundance']
        
    return (abd, metadata)

def plot_pathogens_abundance(pathogens, h5_file, fig_dir=None, thresh=100, factor_name=None):
    '''
    1st col: OTU
    2nd col: abundance
    3rd col: hit on pathogenic db
    4-last col: taxonomy
    '''

    (shared, taxonomy, _, metadata) = load_h5(h5_file)

    pathogenic_otus = np.intersect1d(pathogens.index, taxonomy.index)

    data = get_clustered_zscores(shared[pathogenic_otus].copy(),
                                 factor=factor_name,
                                 metadata=metadata)

    (abd, metadata) = get_shared_info(shared[pathogenic_otus],
                                      factor=factor_name,
                                      metadata=metadata)

    if not factor_name:
        factor_name = 'sampleID'

    suffix = f'_by-{factor_name}'

    # data = data.merge(abd).set_index(factor_name)
    data = data.set_index(factor_name)

    cols = np.intersect1d(metadata.columns,
                          ['n_samples', 'Site', 'Month', 'Year', 'Flowpath', 'Aquifer', 'Divide'])

    info = {
        'otu': taxonomy.loc[pathogenic_otus].assign(hit=pathogens.loc[pathogenic_otus]),
        'sample': metadata.loc[np.isin(metadata.index, data.index.unique()), cols]
    }

    clustermap(data, info, fig_dir=fig_dir, thresh=thresh,
               title='biclustered_heatmap_pathogens{}'.format(suffix))

if __name__ == '__main__':
    cfg = parse_config()
    args = parse_args()

    sp_file = Path(cfg.get('microbial', 'species'))
    db_file = Path(cfg.get('misc', 'patho_db'))
    h5_file = cfg.get('misc', 'h5')
    fig_dir = cfg.get('misc', 'fig_dir')
    thresh = cfg.get('misc', 'otu_thresh')

    pathogens = get_pathogens(sp_file, db_file)

    plot_pathogens_abundance(pathogens, h5_file,
                             fig_dir=fig_dir, thresh=thresh, factor_name=args.factor)
