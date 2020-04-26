from bokeh.plotting import figure

from ecotools.api.bokeh_util import bokeh_legend_out, bokeh_save, bokeh_facets, get_palette
from ecotools.api.bokeh_util import TOOLS, PADDING

@bokeh_save
@bokeh_facets
@bokeh_legend_out
def taxa_barplot(metagenome, output=None, factor=None, cmap={},
                 taxa_file=None, taxa=None, clade=False,
                 norm=False, rank=None, thresh=None, top=-1):

    metagenome = metagenome.copy()
    metagenome.preprocess(factor=factor, norm=norm, rank=rank, top=top,
                          taxa_file=taxa_file, taxa=taxa, clade=clade,)

    # Discard OTU with a min relative abundance below 1%
    main_otus = metagenome.abundance.get_most_abundant_otus(thresh=thresh)

    data = metagenome.abundance.data[main_otus]

    data.columns = (data.columns
                    .str.replace('[/\-]', '_')
                    .str.replace('(\(.*\))', '')
                    .str.slice(stop=100))

    tooltips = zip(data.columns, '@'+data.columns+'{0.00%}'*norm)

    factors = data.index.tolist()
    # if all(f.isdigit() for f in factors):
    #     factors = sorted(factors, key=int)

    width = max(600, metagenome.n_samples()*50) + PADDING
    height = 400 + PADDING

    taxa_to_add = set(data.columns).difference(set(cmap))
    add_cmap = zip(list(taxa_to_add), get_palette(len(taxa_to_add)))
    cmap.update(dict(add_cmap))

    colors = [cmap[x] for x in data.columns]

    p = figure(plot_height=height, width=width, min_border=PADDING, x_range=factors,
               tools=TOOLS, tooltips=list(tooltips),
               title=output.stem)

    p.vbar_stack(data.columns, x=data.index.name, source=data.reset_index(),
                 width=.8, color=colors,
                 legend_label=data.columns.tolist())

    p.xaxis.major_label_orientation = 'vertical'
    
    return p
