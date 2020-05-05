from pathlib import Path

import pandas as pd

from ecotools.rpy2_util import metamds
from ecotools.util import guess_subsampling_level
from ecotools.plotting.scatter import scatterplot

def nmds(metagenome, metric='braycurtis', k=3, trymax=200, parallel=3, subsample=True,
         cache=False, inplace=False,
         plot=False, **plot_kw):

    nmds_file = Path(metagenome.outdir, 'NMDS_{}.csv'.format(metric))

    if not inplace:
        metagenome = metagenome.copy()

    if metagenome.distance_matrix is None:
        if subsample:
            subsampling_level = guess_subsampling_level(metagenome.sample_sizes())
            metagenome.subsample(level=subsampling_level)

        if 'unifrac' in metric.lower():
            metagenome.seq_data.update_tree(metagenome.otus(), app='FastTree')
            metagenome.seq_data.update_tree(metagenome.otus(), app='FastTree')
        else:
            metagenome.wisconsin()
            tree = None

        metagenome.compute_distance_matrix(metric=metric, tree=tree, cache=cache)

    if nmds_file.is_file() and cache:
        nmds_components = pd.read_csv(nmds_file, index_col=0)
    else:
        nmds_components = metamds(metagenome.distance_matrix)
        nmds_components.to_csv(nmds_file)
    
    if plot:
        plot_kw['output'] = 'NMDS_{}_by_{}.html'.format(metric, plot_kw['hue'])
        scatterplot(metagenome, nmds_components, **plot_kw)
    
    return nmds_components

