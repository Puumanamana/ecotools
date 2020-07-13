from itertools import combinations
import pandas as pd
import numpy as np

import rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

rpy2.rinterface_lib.callbacks.consolewrite_warnerror = lambda x: None

def pandas_to_r(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_from_pd_df = ro.conversion.py2rpy(df)

    return r_from_pd_df

def r_to_pandas(df):
    with localconverter(ro.default_converter + pandas2ri.converter):
        df_py = ro.conversion.rpy2py(df)

    return df_py

def pandas_to_distR(dist_df):
    stats_pkg = importr('stats')
    dists_r = stats_pkg.dist(pandas_to_r(dist_df))

    return dists_r


def to_phyloseq(abundance, taxonomy, metadata, tree_path=None):
    phylo = importr('phyloseq')
    base = importr('base')

    args = [
        phylo.otu_table(pandas_to_r(abundance), taxa_are_rows=False),
        phylo.tax_table(base.as_matrix(pandas_to_r(taxonomy))),
        phylo.sample_data(pandas_to_r(metadata))
    ]
    if tree_path is not None:
        ape = importr('ape')
        args.append(ape.read_tree(tree_path))

    phylo_obj = phylo.phyloseq(*args)

    return phylo_obj

def phyloseq_ordinate(obj, method='PCoA', distance='bray', covariates=None, k=2):
    if method.lower() == 'pcoa':
        method = 'PCoA'
    else:
        method = method.upper()
    
    phylo = importr('phyloseq')

    kw = dict(method=method, distance=distance, trymax=500)
    if method == 'CCA':
        kw['formula'] = ro.Formula('. ~ {}'.format('+'.join(covariates)))

    model = phylo.ordinate(obj, **kw)

    if method == 'PCoA':
        ordination = dict(
            inertia=['{:.2%}'.format(x) for x in model.rx2('values').rx2('Relative_eig')[:k]],
            sample=np.array(model.rx2('vectors'))[:, :k]
        )
    elif method == 'NMDS':
        info = 'stress={:.2f}-converged={}'.format(model.rx2('stress')[0], int(model.rx2('converged')[0]))
        ordination = dict(
            qc=np.concatenate([model.rx2('stress'), model.rx2('converged')]),
            inertia=[info] * k,
            sample=np.array(model.rx2('points')),
            species=np.array(model.rx2('species'))
        )
    elif method == 'CCA':
        vegan = importr('vegan')
        scores = vegan.scores(model, display=ro.StrVector(['sp', 'wa', 'bp']))
        
        ordination = dict(
            inertia=['{:.2%}'.format(x) for x in model.rx2('CCA').rx2('eig')[:k]],
            sample=np.array(scores.rx2('sites')),
            feature=np.array(scores.rx2('species')),
            biplot=np.array(scores.rx2('biplot'))
        )

    return ordination        
    

def vegdist(abundance, metric='bray', binary=False, r_obj=False):
    if metric.lower().startswith('bray'):
        metric = 'bray'
    vegan_pkg = importr('vegan')
    abundance_r = pandas_to_r(abundance.sort_index())
    
    dists_r = vegan_pkg.vegdist(abundance_r, method=metric,
                                binary=binary, diag=False)

    if r_obj:
        return dists_r

    indices = combinations(abundance.index, 2)
    dists = pd.Series(r_to_pandas(dists_r), index=pd.Index(indices))
    diag = pd.Series(1, index=pd.Index(zip(abundance.index, abundance.index)))
    
    dists = pd.concat([dists, diag]).sort_index().rename(metric)
    return dists

def metamds(dists_r, metric='bray', k=2, trymax=500, parallel=20):
    vegan_pkg = importr('vegan')
    nmds_obj = vegan_pkg.metaMDS(dists_r, distance=metric, k=k,
                                 trymax=trymax, parallel=parallel)
    nmds_components = vegan_pkg.scores(nmds_obj)

    nmds_components = pd.DataFrame(
        np.array(nmds_components),
        columns=["nmds_{}".format(i+1) for i in range(k)],
    )

    stress = np.array(nmds_obj.rx2('stress'))[0]

    import ipdb;ipdb.set_trace()
    return (stress, nmds_components)

def pcoa_ape(dists_r, k=2):

    ape_pkg = importr('ape')
    pcoa_obj = ape_pkg.pcoa(dists_r)

    components = pd.DataFrame(
        r_to_pandas(pcoa_obj.rx2('vectors'))[:, :k],
        columns=["pcoa_{}".format(i+1) for i in range(k)],
    )

    return components

def pcoa_phylo(phylo_obj, k=2):

    phylo_pkg = importr('phyloseq')
    pcoa_obj = phylo_pkg.ordinate(phylo_obj, method='PCoA')

    components = pd.DataFrame(
        np.array(pcoa_obj.rx2('vectors'))[:, :k],
        columns=["pcoa_{}".format(i+1) for i in range(k)],
    )

    explained = pd.Series(
        np.array(pcoa_obj.rx2('values').rx2('Relative_eig'))[:k],
        index=["pcoa_{}".format(i+1) for i in range(k)]
    )
    
    return (explained, components)


def cca_vegan(abundance, metadata, k=2):
    (x, y) = (pandas_to_r(abundance), pandas_to_r(metadata))
    
    vegan_pkg = importr('vegan')
    model = vegan_pkg.cca(x, y, k=k)
    scores = vegan_pkg.scores(model, display=ro.StrVector(['sp', 'wa', 'bp']))

    constr_inertia = (np.array(model.rx2('CCA').rx2('tot.chi')) /
                      np.array(model.rx2('tot.chi')))[0]

    labels = ['cca{}'.format(i) for i in range(2)]

    components = {
        'biplot': pd.DataFrame(
        np.array(scores.rx2('biplot')),
        index=metadata.columns, columns=labels
        ),
        'sample': pd.DataFrame(
            np.array(scores.rx2('sites')),
            index=abundance.index, columns=labels
        ),
        'feature': pd.DataFrame(
            np.array(scores.rx2('species')),
            index=abundance.columns, columns=labels
        )}

    return (constr_inertia, components)
    

def permanova_r(dists_r, factors, permutations=9999):
    factors_r = pandas_to_r(factors)
     
    pwAdonis_pkg = importr('pairwiseAdonis')

    r_model = pwAdonis_pkg.pairwise_adonis(dists_r, factors_r,
                                           perm=permutations)

    r_res_items = list(r_model.items())
    index = pd.Index(r_res_items[0][1].levels, name=r_res_items[0][0])

    model_results = pd.DataFrame(
        dict(r_res_items[1:]),
        columns=r_model.colnames[1:],
        index=index
    ).rename(columns={'p.adjusted': 'pval_adj', 'F.Model': 'statistic'})     

    return model_results
