from bokeh.plotting import figure

from ecotools.api.util import bokeh_legend_out, bokeh_save, get_palette
from ecotools.api.util import TOOLS, PADDING

@bokeh_save
@bokeh_legend_out
def taxa_barplot(metagenome, output=None, factor=None,
                 taxa_file=None, taxa=None, clade=False,
                 norm=False, rank=None, top=-1, padding=200):

    metagenome = metagenome.copy()
    metagenome.preprocess(factor=factor, norm=norm, rank=rank, top=top,
                          taxa_file=taxa_file, taxa=taxa, clade=clade,)

    # Discard OTU with a min relative abundance below 0.1%
    main_otus = metagenome.abundance.get_most_abundant_otus(thresh=0.001)
    data = metagenome.abundance.data[main_otus]

    data.columns = (data.columns
                    .str.replace('[/\-]', '_')
                    .str.replace('(\(.*\))', '')
                    .str.slice(stop=100))

    tooltips = zip(data.columns, '@'+data.columns+'{0.00%}'*norm)

    factors = data.index.tolist()
    if all(f.isdigit() for f in factors):
        factors = sorted(factors, key=int)

    width = max(800, metagenome.n_samples()*100) + PADDING
    height = 600 + PADDING

    palette = get_palette(data.shape[1])

    p = figure(plot_height=height, width=width, min_border=padding,
               x_range=factors,
               tools=TOOLS, tooltips=list(tooltips),
               title=output.stem)

    p.vbar_stack(data.columns, x=data.index.name, source=data.reset_index(),
                 width=.6, color=palette[:metagenome.n_otus()],
                 legend_label=data.columns.tolist())

    p.xaxis.major_label_orientation = 'vertical'
    
    return p
