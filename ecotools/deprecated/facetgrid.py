import numpy as np
import pandas as pd

from bokeh.layouts import gridplot
from bokeh.io import output_file, save
from bokeh.models import Legend, Title, HoverTool

from ecotools.util import get_palette
from ecotools.parsing import parse_config

CFG = parse_config()['bokeh']

def format_legend(plot, leg_order=None, shorten=True):

    leg_it = plot.legend.items[::-1]
    
    txt_len = CFG['leg_txt_len']

    if shorten and 'value' in leg_it[0].label:
        for i, item in enumerate(leg_it):
            lab_i = item.label['value']
            leg_it[i].label['value'] = (lab_i[:txt_len] + '..'
                                        if len(lab_i) > txt_len
                                        else lab_i)
            
    legends = []

    col = 0
    while leg_it:
        if col == CFG['leg_maxcol']:
            items = [leg_it.pop() for _ in range(len(leg_it))]
        else:
            items = [leg_it.pop() for _ in range(CFG['leg_nrows']) if leg_it]

        legend = Legend(items=items, location=(20, 0), glyph_height=20, glyph_width=20,
                        padding=0, spacing=0, border_line_color=None,
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
                 row_order=None, col_order=None, hue_order=None, paired_colors=False,
                 width=600, height=600, scale=1,
                 palette='Turbo256', randomize_palette=False, outdir='.'):
        self.data = data.copy()
        self.col = col
        self.row = row
        self.hue = hue
        self.ordering = {'row': row_order, 'col': col_order, 'hue': hue_order}
        self.height = int(scale*height)
        self.width= int(scale*width)
        self.ncols = 1
        self.nrows = 1
        self.cmap = None
        self.paired_colors = paired_colors
        self.palette = palette
        self.randomize_palette = randomize_palette
        self.outdir = outdir
        self.plots = []

    def update_cmap(self, data):
        if self.hue is None:
            return

        if self.ordering.get('hue', None) is None:
            self.ordering['hue'] = sorted(data[self.hue].unique())
    
        self.data[self.hue] = pd.Categorical(self.data[self.hue],
                                             self.ordering['hue'])
        self.data.sort_values(by=self.hue, inplace=True)

        if self.cmap is None:
            colors = get_palette(len(self.ordering['hue']), self.palette,
                                 random=self.randomize_palette, paired=self.paired_colors)
            self.cmap = dict(zip(self.ordering['hue'], colors))

            for label in self.cmap:
                if label.startswith('Others'):
                    self.cmap[label] = '#cccccc'

    def get_row_col_combs(self):
        factors = {}

        if self.row is None and self.col is None:
            return ([True], np.ones(len(self.data), dtype=bool))
        
        if self.row is not None:
            factors[self.row] = self.data[self.row]
            self.nrows = len(factors[self.row].unique())

            if self.ordering['row'] is None:
                self.ordering['row'] = sorted(factors[self.row].unique())
                
        if self.col is not None:
            factors[self.col] = self.data[self.col]
            self.ncols = len(factors[self.col].unique())

            if self.ordering['col'] is None:
                self.ordering['col'] = sorted(factors[self.col].unique())

        factor_values = pd.Series(zip(*[fact for fact in factors.values()]),
                                  index=self.data.index)
        factor_combs = pd.MultiIndex.from_product(
            [x for k, x in self.ordering.items() if x is not None and k != 'hue'],
            names=list(factors.keys()))

        return (factor_combs, factor_values)

    def map(self, func, x=None, y=None, tooltips=[], rotate=0, **kwargs):
        is_new = (len(self.plots) == 0)

        if not isinstance(x, list):
            x = [x]

        if self.hue is not None:
            x += [self.hue] 
            self.update_cmap(self.data)
            self.data['color'] = [self.cmap[lvl] for lvl in self.data[self.hue]]
            kwargs['fill_color'] = 'color'
            kwargs['hue_order'] = self.ordering['hue']

        (factor_combs, factor_values) = self.get_row_col_combs()        
        for i, pair in enumerate(factor_combs):
            data = self.data.loc[(factor_values == pair)].copy()
            
            if data.size == 0:
                self.plots.append(None)
                continue

            prev_plot = None
            if not is_new: prev_plot = self.plots[i]
            
            if self.hue is not None:
                tooltips = np.append(tooltips, self.hue)
                if is_new:
                    kwargs['legend_field'] = self.hue

            x = [xi for xi in x if xi is not None]
            if len(x) == 1 and 'legend_field' in kwargs:
                kwargs.pop('legend_field')
                
            p = func(x=x, y=y, data=data, p=prev_plot,
                     width=self.width, height=self.height,
                     name=func.__name__, **kwargs)

            if rotate != 0:
                p.xaxis.major_label_orientation = rotate

            tooltips_fmt = pd.Series(tooltips, dtype='object').drop_duplicates()
            tooltips_fmt = list(zip(tooltips_fmt, '@'+tooltips_fmt))
            hover_tool = HoverTool(tooltips=tooltips_fmt, names=[func.__name__])
            hover_tool.point_policy='snap_to_data'
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

    def simple_map(self, fn, **kwargs):
        for i, p in enumerate(self.plots):
            p = fn(p=p, **kwargs)
            self.plots[i] = p

    def save(self, filename):

        first_p = next(p for p in self.plots if p is not None)
            
        grid = gridplot(self.plots, ncols=self.ncols,
                        plot_width=first_p.plot_width,
                        plot_height=first_p.plot_height)

        output_file('{}/{}'.format(self.outdir, filename))
        save(grid)
        
