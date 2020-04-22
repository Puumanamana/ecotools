from bokeh.plotting import figure

from util import bokeh_legend_out, bokeh_save, get_palette, TOOLS

@bokeh_save
@bokeh_legend_out
def taxa_barplot(metagenome, output=None, factor=None, taxa_file=None, norm=False, rank=None, top=-1, padding=200):

    metagenome = metagenome.copy()
    metagenome.preprocess(factor=factor, taxa_file=taxa_file, norm=norm, rank=rank, top=top)

    padding = 200
    width = padding + max(800, metagenome.n_samples()*100)
    height = 500 + padding

    palette = get_palette(metagenome.n_otus())

    data = metagenome.abundance.data
    data.columns = data.columns.str.replace('/', '_').str.slice(stop=100)

    tooltips = zip(data.columns, '@'+data.columns+'{0.00%}'*norm)

    factors = data.index.tolist()
    if all(f.isdigit() for f in factors):
        factors = sorted(factors, key=int)

    p = figure(plot_height=height, width=width, min_border=padding,
               x_range=factors,
               tools=TOOLS, tooltips=list(tooltips),
               title=output.stem)

    p.vbar_stack(data.columns, x=data.index.name, source=data.reset_index(),
                 width=.3, color=palette[:metagenome.n_otus()],
                 legend_label=data.columns.tolist())

    p.xaxis.major_label_orientation = 'vertical'
    
    return p
