import textwrap
import time
from itertools import combinations

import pandas as pd

from ecotools.parsing import parse_config

CFG = parse_config()['bokeh']

def timer(func):
    def wrapper(*args, **kw):
        start = time.time()
        result = func(*args, **kw)
        end = time.time()

        print(textwrap.dedent(f'''
        ========{func.__name__}========
        Elapsed time: {end-start:.5f}s'''))
        
        return result
    return wrapper

def pairwise(func):
    def wrapper(*args, **kwargs):

        (metagenome, factor_name) = args[:2]
        
        factor = metagenome.metadata.factor_data(factor_name)

        if 'pairwise' not in kwargs or not kwargs['pairwise']:
            result = func(metagenome, factor, **kwargs)
            return pd.concat(result)

        results = {}
        for (lvl1, lvl2) in combinations(factor.unique(), 2):
            mg = metagenome.copy()
            subset = mg.abundance.index[(factor == lvl1) | (factor == lvl2)]
            mg.subset_samples(subset)

            name = '{}-vs-{}'.format(*sorted([lvl1, lvl2]))
            results[name] = func(mg, factor[subset], *args[2:], **kwargs)

        return pd.concat(results).reset_index()

    return wrapper


def strata(func):
    def wrapper(*args, **kwargs):
        metagenome = args[0]

        strata_names = kwargs.pop('strata', None)

        if not isinstance(strata_names, list):
            strata_names = [strata_names]
        
        strata_names = [x for x in strata_names if x is not None]
        
        if not strata_names:
            return func(*args, **kwargs)

        sample_strata = metagenome.metadata.factor_data(strata_names).apply(tuple, axis=1)
        strata_combs = sample_strata.drop_duplicates().tolist()

        results = {}
        for strata in strata_combs:
            subset = metagenome.abundance.index[sample_strata == strata]
            level_name = '-'.join(strata)

            if subset.size > 1:            
                mg = metagenome.copy()
                mg.subset_samples(subset)

                kwargs['strata'] = level_name
                result = func(mg, *args[1:], **kwargs)
                results[level_name] =  result

        return results

    return wrapper
