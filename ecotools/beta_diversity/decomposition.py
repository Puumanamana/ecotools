from pathlib import Path

import pandas as pd

from ecotools.decorators import strata
from ecotools.rpy2_util import metamds


@strata
def nmds(metagenome, metric='braycurtis', k=3, trymax=500, parallel=3,
         subsample=True, subsampling_level=-1,
         cache=False, inplace=False, **kw):

    suffix = '_'.join([metric, kw.pop('strata', [''])[0]]).strip('_')
    nmds_file = Path(metagenome.outdir, 'NMDS_{}.csv'.format(suffix))

    if nmds_file.is_file() and cache:
        nmds_components = pd.read_csv(nmds_file, index_col=0)
        return nmds_components

    if metagenome.distance_matrix is None:
        print('Could not find the distance matrix. Did you compute it first?')
        return
    
    if not inplace:
        metagenome = metagenome.copy()

    distance_matrix = metagenome.distance_matrix.loc[metagenome.abundance.index,
                                                     metagenome.abundance.index]

    nmds_components = metamds(distance_matrix, trymax=trymax, k=k, parallel=parallel)
    nmds_components.to_csv(nmds_file)
    
    return nmds_components    
