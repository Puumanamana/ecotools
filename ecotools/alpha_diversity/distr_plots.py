from ecotools.plotting.barplot import stacked_barplot
from ecotools.plotting.boxplot import boxplot

def taxa_barplot(metagenome, taxa_files=None, taxa=None, clade=False,
                 relabund=False, rank=None, thresh=None, top=-1, **plot_kw):

    metagenome = metagenome.copy()
    metagenome.preprocess(relabund=relabund, rank=rank, top=top,
                          taxa_files=taxa_files, taxa=taxa, clade=clade)

    # Discard OTU with a min relative abundance below 1%
    main_otus = metagenome.abundance.get_most_abundant_otus(thresh=thresh)

    metagenome.subset_otus(otus=main_otus)

    metagenome.abundance.data.columns = metagenome.abundance.data.columns.str.replace('[\(\)\-]', '_')
    metagenome.taxonomy.data.index = metagenome.taxonomy.data.index.str.replace('[\(\)\-]', '_')
    stacked_barplot(metagenome, randomize_palette=True, **plot_kw)

    
def richness_distr(metagenome, metric='n_otus', x=None, hue=None,
                   taxa_files=None, taxa=None, clade=False, relabund=True, rank=None, **plot_kw):

    metagenome = metagenome.copy()
    metagenome.preprocess(relabund=relabund, rank=rank,
                          taxa_files=taxa_files, taxa=taxa, clade=clade)
    metagenome.compute_alpha_diversity(metric)

    boxplot(metagenome, x=x, hue=hue, **plot_kw)
    

    

    

