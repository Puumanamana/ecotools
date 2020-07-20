from pathlib import Path

import pandas as pd

from ecotools.decorators import strata, timer
from ecotools.rpy2_util import phyloseq_ordinate 

@strata
@timer
def ordinate(mg, subsample=False, otu_subset=None, covariates=None, k=2,
             cache=False, strata='', sep=','):
    
    if subsample:
        mg.subsample(plot=True)

    output = mg.io['ordination']((mg.run['clade'], strata))

    if cache and output.is_file():
        components = pd.read_csv(output, dtype={0: str})
        components = components.set_index(components.columns[0])
        return components

    if otu_subset is not None:
        if isinstance(otu_subset, Path):
            mg.subset_otus(otu_subset, sep=sep)

    phylo_obj = mg.to_phyloseq()
    result = phyloseq_ordinate(phylo_obj, mg.run['ordination'], mg.run['metric'],
                               covariates=covariates)

    sample_cols = []
    for i in range(k):
        label = '{}_{}'.format(mg.run['ordination'], i)
        # if 'inertia' in result:
        #     label += ' ({})'.format(result['inertia'][i])
        sample_cols.append(label)
    
    result['sample'] = pd.DataFrame(result['sample'], index=mg.index, columns=sample_cols)

    for key in ['sample', 'feature', 'biplot']:
        if strata and key in result:
            result[key]['strata'] = strata

    if 'feature' in result:
        result['feature'] = pd.DataFrame(
            result['feature'], index=mg.columns,
            columns=['{}_{}'.format(mg.run['ordination'], i) for i in range(k)])

    if 'biplot' in result:
        result['biplot'] = pd.DataFrame(
            result['biplot'], index=covariates,
            columns=['{}_{}'.format(mg.run['ordination'], i) for i in range(k)])

    return result
