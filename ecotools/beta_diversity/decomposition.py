import sys
from pathlib import Path

import pandas as pd

from ecotools.decorators import strata, timer
from ecotools.rpy2_util import metamds, pcoa_ape, pandas_to_distR


@strata
@timer
def decompose(mg, decomposition='nmds', metric='bray',
              cache=False, strata='', return_distances=False,
              **decomposition_kw):
    
    output = Path(mg.outdir, f'{decomposition}_{metric}_{strata}.csv')

    if cache and output.is_file():
        nmds_components = pd.read_csv(output, dtype={0: str})
        return nmds_components.set_index(nmds_components.columns[0])

    if mg.distance_matrix is None:
        mg.compute_distance_matrix(metric=metric, vegan=True)

    dists_r = pandas_to_distR(mg.distance_matrix.unstack())
    if decomposition.lower() == 'nmds':
        components = metamds(dists_r, **decomposition_kw)
    elif decomposition.lower() == 'pcoa':
        components = pcoa_ape(dists_r)
    else:
        sys.exit(f'Unknown decomposition method {decomposition}')
        
    components.index = mg.index

    # Save results:
    components.to_csv(output)
    
    return components    
