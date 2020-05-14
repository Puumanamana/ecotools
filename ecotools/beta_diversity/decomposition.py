from pathlib import Path

import pandas as pd

from ecotools.decorators import strata
from ecotools.rpy2_util import metamds
from ecotools.util import guess_subsampling_level

@strata
def nmds(metagenome, metric='braycurtis', k=3, trymax=500, parallel=3, subsample=True,
         cache=False, inplace=False, **kw):

    suffix = '_'.join([metric, kw.pop('strata', [''])[0]]).strip('_')
    nmds_file = Path(metagenome.outdir, 'NMDS_{}.csv'.format(suffix))

    if nmds_file.is_file() and cache:
        nmds_components = pd.read_csv(nmds_file, index_col=0)
        return nmds_components

    if not inplace:
        metagenome = metagenome.copy()

    if metagenome.distance_matrix is None:
        if subsample:
            subsampling_level = guess_subsampling_level(metagenome.sample_sizes())
            metagenome.subsample(level=subsampling_level)

        if 'unifrac' in metric.lower():
            metagenome.seq_data.update_tree(metagenome.otus(), app='FastTree')
            metagenome.seq_data.update_tree(metagenome.otus(), app='FastTree')
            tree = metagenome.seq_data.load_tree()
        else:
            metagenome.normalize('wisconsin')
            tree = None

        metagenome.compute_distance_matrix(metric=metric, tree=tree, cache=cache)

    distance_matrix = metagenome.distance_matrix.loc[metagenome.samples(), metagenome.samples()]

    nmds_components = metamds(distance_matrix, trymax=trymax, k=k, parallel=parallel)
    nmds_components.to_csv(nmds_file)
    
    return nmds_components    
