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

    def __init__(self, data=None,
                 hue=None, col=None, row=None,
                 hue_order=None, col_order=None, row_order=None,
                 width=600, height=600, scale=1, col_wrap=None,
                 palette='Turbo256', randomize_palette=False, paired_colors=False, outdir='.'):
        self.data = data.copy()
        self.cols = col
        self.rows = row
        self.hue = hue
        self.col_wrap = col_wrap
        self.height = 'auto' if height=='auto' else int(scale*height)
        self.width = 'auto' if width=='auto' else int(scale*width)
        self.cmap = None
        self.colors = {'palette': palette, 'random': randomize_palette, 'paired': paired_colors}
        self.outdir = outdir
        self.plots = []
        
        self.format_data(rows=row_order, cols=col_order, hue=hue_order)
        self.update_cmap()

    def format_data(self, **kwargs):
        '''
        Add fake columns in case no hue, col or row arguments is supplied
        '''
        for (attr, categories) in kwargs.items():
            if getattr(self, attr) is None:
                setattr(self, attr, ' ')
                self.data[' '] = ' '

            name = getattr(self, attr)
            if categories is None:
                categories = self.data[name].sort_values().unique()
                            
            self.data[name] = pd.Categorical(self.data[name], categories)

        self.data.sort_values(by=self.hue, inplace=True)
        
    def update_cmap(self):
        if self.cmap is None:
            categories = self.data[self.hue].cat.categories
            colors = get_palette(len(categories), **self.colors)
            self.cmap = dict(zip(categories, colors))

            for label in self.cmap:
                if label.startswith('Others'):
                    self.cmap[label] = '#cccccc'

    def get_row_col_combs(self):
        '''
        Get boolean vectors to subset rows and columns
        ''' 
        factor_values = self.data[[self.rows, self.cols]].apply(tuple, axis=1)
        factor_combs = pd.MultiIndex.from_product(
            [self.data[self.rows].cat.categories, self.data[self.cols].cat.categories],
            names=[self.rows, self.cols]
        )

        return (factor_combs, factor_values)
        
    def map(self, func, x=None, y=None, tooltips=[], rotate=0, **kwargs):
        is_first_map = (len(self.plots) == 0)

        if not isinstance(x, list):
            x = [x]

        if self.hue != ' ':
            x += [self.hue] 
            self.data['color'] = [self.cmap[lvl] for lvl in self.data[self.hue]]
            kwargs['fill_color'] = 'color'
            kwargs['hue_order'] = self.data[self.hue].cat.categories

        (factor_combs, factor_values) = self.get_row_col_combs()        
        for i, pair in enumerate(factor_combs):
            data = self.data.loc[(factor_values == pair)].copy()
            
            if data.size < CFG['min_per_facet']:
                self.plots.append(None)
                continue

            prev_plot = None
            if not is_first_map: prev_plot = self.plots[i]
            
            if self.hue != ' ':
                # Add legend
                tooltips = np.append(tooltips, self.hue)
                if is_first_map:
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
            if is_first_map and pair != (' ', ' '):
                subtitle = Title(
                    text=', '.join('{}: {}'.format(*it) for it in zip(factor_combs.names, pair)
                                   if it[1] != ' '),
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

        ncols = len(self.data[self.cols].cat.categories)
        if self.col_wrap is not None:
            ncols = self.col_wrap
            
        grid = gridplot(self.plots, ncols=ncols,
                        plot_width=first_p.plot_width, plot_height=first_p.plot_height)

        output_file('{}/{}'.format(self.outdir, filename))
        save(grid)
        
