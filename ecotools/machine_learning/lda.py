import pickle
from sklearn.decomposition import LatentDirichletAllocation
import pandas as pd

import pyLDAvis

from ecotools.plotting.heatmap import clustermap
from ecotools.decorators import strata

def ldavis_show(metagenome, sample_probs, otu_probs, output=None):

    taxa_info = (metagenome.taxonomy
                 .data.loc[metagenome.otus(), ['Class', 'Genus']]
                 .apply(lambda x: '_'.join(x), axis=1))

    LDAvis_prepared = pyLDAvis.prepare(
        otu_probs, # (topics x otus)
        sample_probs, # (samples x topics)
        metagenome.abundance.data.sum(axis=1), # (samples)
        taxa_info, # (otus)
        metagenome.abundance.data.sum(axis=0)) # (samples)

    LDAvis_data_filepath = '{}/ldavis_prep.pkl'.format(str(metagenome.outdir))

    with open(LDAvis_data_filepath, 'wb') as f:
        pickle.dump(LDAvis_prepared, f)
        
    # load the pre-prepared pyLDAvis data from disk
    with open(LDAvis_data_filepath, 'rb') as f:
        LDAvis_prepared = pickle.load(f)
    pyLDAvis.save_html(LDAvis_prepared, '{}/{}'.format(metagenome.figdir, output))
    
def sample_topics_clustermap(metagenome, topics, **kwargs):
    metadata = metagenome.metadata.factor_data()
    data = pd.concat([topics, metadata], axis=1).rename_axis(index='groups').reset_index()
    data = data.melt(id_vars=list(metagenome.metadata.qual_vars) + ['groups'],
                     var_name='topics', value_name='weight')

    data.set_index(['groups', 'topics'], inplace=True)

    clustermap(metagenome, data, value_col='weight', standardize=False, **kwargs)

def otu_topics_clustermap(metagenome, topics, **kwargs):
    top_otus = topics.index[topics.max(axis=1) > 0.01]
    
    tax = metagenome.taxonomy.data.drop(columns='Species')
    data = (pd.concat([topics, tax], axis=1).loc[top_otus]
            .rename_axis(index='otus')
            .reset_index())
    data = data.melt(id_vars=list(tax.columns) + ['otus'],
                     var_name='topics', value_name='weight')

    data.set_index(['otus', 'topics'], inplace=True)

    clustermap(metagenome, table=data, value_col='weight', standardize=False, **kwargs)

@strata
def lda_model(metagenome, groups=None, k=5, plot=False, strata=None, **kwargs):

    metagenome = metagenome.copy()
    metagenome.taxonomy.clean_labels()

    if groups is not None:
        metagenome.group_samples(groups)

    model = LatentDirichletAllocation(n_components=k)
    model.fit(metagenome.abundance.data)

    sample_topics = pd.DataFrame(
        model.transform(metagenome.abundance.data),
        index=metagenome.samples(),
        columns=['topic_{}'.format(i+1) for i in range(k)]
    )

    otu_topics = pd.DataFrame(
        model.components_ / model.components_.sum(axis=1)[:, None],
        index=sample_topics.columns,
        columns=metagenome.otus()
    )

    if plot:
        suffix = ''
        if strata is not None:
            suffix = '_{}'.format('-'.join(strata))
            
        sample_topics_clustermap(metagenome, sample_topics,
                                 output='lda_sample-topics{}.html'.format(suffix), **kwargs)
        otu_topics_clustermap(metagenome, otu_topics.T, dim='feature',
                              output='lda_otu-topics{}.html'.format(suffix), **kwargs)        
        # ldavis_show(metagenome, sample_topics, otu_topics,
        #             output='ldavis{}.html'.format(suffix))
    


    
