import textwrap

import numpy as np
import pandas as pd

from bokeh.layouts import gridplot
from bokeh.io import output_file, save
from bokeh.models import Legend, Title, HoverTool

from ecotools.util import get_palette
from ecotools.parsing import parse_config

CFG = parse_config()['bokeh']

def format_legend(plot, shorten=True):

    leg_items = plot.legend.items[::-1]

    if shorten and 'value' in leg_items[0].label:
        for i, item in enumerate(leg_items):
            leg_items[i].label['value'] = textwrap.shorten(
                item.label['value'], CFG['leg_txt_len']
            )
    legends = []

    col = 0
    while leg_items:
        if col == CFG['leg_maxcol']:
            items = [leg_items.pop() for _ in range(len(leg_items))]
        else:
            items = [leg_items.pop() for _ in range(CFG['leg_nrows']) if leg_items]

        legend = Legend(items=items, location=(25, 0), glyph_height=25, glyph_width=25,
                        padding=0, spacing=0,
                        border_line_color=None,
                        label_text_font_size=CFG['leg_fs'])

        col += 1
        legends.append(legend)

    return legends

def format_tooltips(tips):
    if tips is None:
        return
    if len(tips[0]) == 2:
        return tips

    tips = [(x, '@{}'.format(x)) for x in tips]
    return tips

class BokehFacetGrid:

    def __init__(self, hue=None, col=None, row=None, data=None,
                 width=500, height=500, scale=1, randomize_palette=False, outdir='.'):
        self.data = data
        self.col = col
        self.row = row
        self.hue = hue
        self.height = scale*height
        self.width= scale*width
        self.ncols = 1
        self.nrows = 1
        self.cmap = {}
        self.randomize_palette = randomize_palette
        self.outdir = outdir
        self.plots = []

    def update_cmap(self, data):
        if self.hue is None:
            return
        uniq_hue = data[self.hue].unique()
        new_hue_values = [x for x in uniq_hue if x not in self.cmap]
        colors = get_palette(len(new_hue_values), random=self.randomize_palette)
        self.cmap.update(dict(zip(new_hue_values, colors)))

    def get_row_col_combs(self):
        factors = {}
        if self.row is None and self.col is None:
            return ([True], np.ones(len(self.data), dtype=bool))
        if self.row is not None:
            factors[self.row] = self.data[self.row]
            self.nrows = len(factors[self.row].unique())
        if self.col is not None:
            factors[self.col] = self.data[self.col]
            self.ncols = len(factors[self.col].unique())

        factor_values = pd.DataFrame(factors).apply(tuple, axis=1)
        factor_combs = pd.MultiIndex.from_product([sorted(factors[x].unique()) for x in factors],
                                                  names=list(factors.keys()))
        
        return (factor_combs, factor_values)
        
    def map(self, func, x=None, y=None, tooltips=None, **kwargs):
        (factor_combs, factor_values) = self.get_row_col_combs()
        
        is_new = (len(self.plots) == 0)

        if not isinstance(x, list):
            x = [x]
        if self.hue is not None:
            x += [self.hue]
        
        for i, pair in enumerate(factor_combs):
            data = self.data.loc[(factor_values == pair).tolist()].copy()
            
            if data.size == 0:
                self.plots.append(None)
                continue

            prev_plot = None
            if not is_new:
                prev_plot = self.plots[i]
            
            if self.hue is not None:
                self.update_cmap(data)
                data['color'] = [self.cmap[lvl] for lvl in data[self.hue]]
                kwargs['fill_color'] = 'color'
                data.sort_values(by=self.hue, inplace=True)

                if is_new:
                    kwargs['legend_field'] = self.hue

            p = func(x=x, y=y, data=data, p=prev_plot,
                     width=self.width, height=self.height,
                     name=func.__name__, **kwargs)

            # add hover tips if any
            if tooltips is not None:
                tooltips = pd.Series(tooltips)
                tooltips = list(zip(tooltips, '@'+tooltips))
                hover_tool = HoverTool(tooltips=tooltips, names=[func.__name__])
                p.add_tools(hover_tool)            
            
            # Add facet info in subtitle
            if is_new and pair != True:
                subtitle = Title(
                    text=', '.join('{}: {}'.format(*it) for it in zip(factor_combs.names, pair)),
                    text_font_size=CFG['subtitle_fs'],
                    text_font_style="italic"
                )
                p.add_layout(subtitle, 'above')

            # Format legend and put in outside the plot area
            if 'legend_field' in kwargs:
                legends = format_legend(p)
                p.legend.items.clear()
                for legend in legends:
                    p.add_layout(legend, 'right')

            if len(self.plots) < len(factor_combs):
                self.plots.append(p)

    def save(self, filename):

        grid = gridplot(self.plots, ncols=self.ncols,
                        plot_width=self.plots[0].plot_width,
                        plot_height=self.plots[0].plot_height)

        output_file('{}/{}'.format(self.outdir, filename))
        save(grid)
            
if __name__ == '__main__':

    from ecotools.plotting.scatter import swarmplot
    from ecotools.plotting.barplot import stackplot, barplot
    from ecotools.plotting.boxplot import boxplot

    N = 3000
    yy = np.random.randint(0, 100, N)
    gr = np.random.choice(list("abcdef"), N)
    hu = np.random.choice(list('12'), N)
    co = np.random.choice(['c1', 'c2'], N)

    df = pd.DataFrame(dict(score=yy, group=gr, hue=hu, col=co))

    g = BokehFacetGrid(data=df, width=1200)
    # g.map(stackplot, x=['group', 'col'], y='score')
    g.map(boxplot, x='group', y='score')
    g.map(swarmplot, x='group', y='score', tooltips=['hue', 'group'])
    g.save('test_barplot.html')
