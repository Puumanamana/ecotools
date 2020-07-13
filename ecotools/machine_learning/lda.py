from pathlib import Path
import pickle

from sklearn.decomposition import LatentDirichletAllocation
import pandas as pd

import pyLDAvis

from ecotools.plotting.heatmap import clustermap
from ecotools.plotting.grid import BokehFacetGrid
from ecotools.plotting.scatter import swarmplot
from ecotools.plotting.boxplot import boxplot
from ecotools.decorators import timer

@timer
def lda_model(mg, groups=None, k=5, **kwargs):

    mg = mg.copy()
    mg.taxonomy.clean_labels()

    # if rank is not None:
    #     mg.group_taxa(rank)
    if groups is not None:
        mg.group_samples(groups)

    model = LatentDirichletAllocation(n_components=k, random_state=42)
    model.fit(mg.abundance.data)

    sample_topics = pd.DataFrame(
        model.transform(mg.abundance.data),
        index=mg.abundance.index,
        columns=['topic_{:02}'.format(i+1) for i in range(k)]
    )

    otu_topics = pd.DataFrame(
        (model.components_.T/model.components_.sum(axis=1)).T,
        columns=mg.abundance.columns,
        index=sample_topics.columns
    )

    return {'samples': sample_topics, 'features': otu_topics}


@timer
def lda_boxplot(data, metadata=None, taxonomy=None, x=None, row=None, col=None, rank='Genus', output='lda_plot.html', width=1400, top=10):

    top_otu = (data['features'].stack()
               .rename_axis(index=['topic', 'feature'])
               .rename('weight')
               .reset_index()
               .groupby('topic')
               .apply(lambda x: x.nlargest(top, 'weight').reset_index(drop=True)))

    top_otu[rank] = taxonomy.loc[top_otu.feature, rank].to_numpy()

    top_otu['top_otus'] = (
        top_otu[[rank, 'weight']]
        .apply(lambda x: f'{x.weight:.0%} {x[rank]}', axis=1)
    )

    top_otu = top_otu['top_otus'].unstack()
    top_otu.columns = ['OTU_{}'.format(x+1) for x in top_otu.columns]

    data = pd.concat([data['samples'], metadata], axis=1)
    data = data.melt(id_vars=metadata.columns)
    data = data.merge(top_otu, left_on='variable', right_index=True)

    g = BokehFacetGrid(data=data, hue='variable', row=row, col=col, width=width,
                       outdir=Path(output).parent)
    g.map(boxplot, x=x, y='value', tooltips=top_otu.columns)
    g.map(swarmplot, x=x, y='value', tooltips=metadata.columns)
    g.save(Path(output).name)


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
    
def sample_topics_clustermap(metagenome, topics, output='sample_topics_clustermap.html',
                             row=None, col=None, **kwargs):
    metadata = metagenome.metadata.factor_data()
    data = pd.concat([topics, metadata], axis=1).rename_axis(index='groups').reset_index()
    data = data.melt(id_vars=list(metagenome.metadata.qual_vars) + ['groups'],
                     var_name='topics', value_name='weight')

    g = BokehFacetGrid(data=data, outdir=Path(output).parent, row=row, col=col)
    g.map(clustermap, x='topics', y='groups', z='weight', standardize=False)
    g.save(Path(output).name)

def otu_topics_clustermap(metagenome, topics, output='otu_topics_clustermap.html',
                          col=None, row=None, **kwargs):
    top_otus = topics.index[topics.max(axis=1) > 0.01]
    
    tax = metagenome.taxonomy.data.drop(columns='Species')
    data = (pd.concat([topics, tax], axis=1).loc[top_otus]
            .rename_axis(index='otus')
            .reset_index())

    data = data.melt(id_vars=list(tax.columns) + ['otus'],
                     var_name='topics', value_name='weight')

    g = BokehFacetGrid(data=data, outdir=Path(output).parent, row=row, col=col)
    g.map(clustermap, x='topics', y='otus', z='weight', standardize=False)
    g.save(Path(output).name)
