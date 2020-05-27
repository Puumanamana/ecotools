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
        Ellapsed time: {end-start:.5f}s'''))
        
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

        if 'strata' not in kwargs or kwargs['strata'] is None:
            return func(*args, **kwargs)

        strata_names = kwargs['strata']        
        if not isinstance(strata_names, list):
            strata_names = [strata_names]

        strata_title = '-'.join(strata_names)

        sample_strata = metagenome.metadata.factor_data(strata_names).apply(tuple, axis=1)
        strata_combs = sample_strata.drop_duplicates().tolist()

        results = {}
        for strata in strata_combs:
            subset = metagenome.abundance.index[sample_strata == strata]

            if subset.size > 1:            
                mg = metagenome.copy()
                mg.subset_samples(subset)

                level_name = '-'.join(strata)
                kwargs['strata'] = strata
                result = func(mg, *args[1:], **kwargs)
                results[level_name] =  result

                if isinstance(result, pd.DataFrame):
                    results[level_name][strata_title] = strata[0]
                if is_pd_dict(result):
                    for key in result:
                        results[level_name][key][strata_title] = strata[0]

        if isinstance(result, pd.DataFrame):
            return pd.concat(list(results.values()))
        elif is_pd_dict(result):
            return {k2: pd.concat([
                results[k1][k2] for k1 in results
            ]) for k2 in result}

        return results

    return wrapper

def is_pd_dict(x):
    if isinstance(x, dict):
        return all(isinstance(y, pd.DataFrame) for y in x.values())
    return False

    
         
