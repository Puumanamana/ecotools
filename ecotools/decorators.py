import sys
from pathlib import Path
from itertools import combinations

import pandas as pd

from bokeh.models import Legend
from bokeh.io import output_file, save
from bokeh.layouts import gridplot

from ecotools.parsing import parse_config
from ecotools.util import get_palette

CFG = parse_config()['bokeh']

def bokeh_legend_out(func):
    def wrapper(*args, **kwargs):
        plots = func(*args, **kwargs)
        
        if type(plots) != type([]):
            plots = [plots]

        for p in plots:
            if not p.title.text:
                p.title.text = Path(kwargs['output']).stem.replace('_', ' ')
        
            leg_items = []
            for item in p.legend.items[::-1]:
                if 'value' in item.label:
                    label = item.label['value']
                    item.label['value'] = '{}{}'.format(
                        label[:CFG['leg_txt_len']],
                        '..'*(len(label) > CFG['leg_txt_len']))
                leg_items.append(item)

            p.legend.items.clear()
        
            col = 0
            while leg_items:
                if col == CFG['leg_maxcol']:
                    items = [leg_items.pop() for _ in range(len(leg_items))]
                else:
                    items = [leg_items.pop() for _ in range(CFG['leg_nrows']) if leg_items]
                legend = Legend(items=items, location=(25, 0), glyph_height=25, glyph_width=25)
                p.add_layout(legend, 'right')
                col += 1

                p.legend.label_text_font_size = CFG['leg_fs']
                p.legend.padding = 0
                p.legend.spacing = 0
                p.legend.border_line_color = None
                # p.legend.label_width = LEG_TXT_LEN
        return plots
    return wrapper

def bokeh_save(func):
    def wrapper(*args, **kwargs):
        p = func(*args, **kwargs)

        if 'output' in kwargs and kwargs['output'] is not None:
            output = kwargs['output']
            
            if 'figdir' in args[0].__dict__:
                output = '{}/{}'.format(args[0].__dict__['figdir'], output)
                
            output_file(output)
            save(p)
        else:
            return p

    return wrapper
            

def bokeh_cmap(func):
    '''
    Decorator to merge cmaps between subplots
    '''

    def wrapper(*args, **kwargs):
        
        metagenome = args[0]

        if 'hue' not in kwargs:
            hue_values = args[0].otus()
        else:
            hue_values = metagenome.metadata.data[kwargs['hue']].unique()

        random = False if 'randomize_palette' not in kwargs else kwargs['randomize_palette']
            
        cmap = dict(zip(hue_values, get_palette(len(hue_values), random=random)))
        kwargs['cmap'] = cmap

        return func(*args, **kwargs)

    return wrapper
        

def bokeh_facets(func):
    '''
    Decorator to make facets with factor levels
    '''

    def wrapper(*args, **kwargs):

        if not type(args[0]).__name__ == 'MetagenomicDS':
            sys.exit('First argument must be of type MetagenomicsDS. Found {}'
                     .format(type(args[0]).__name__))

        metagenome = args[0]
        metagenome.sort_otus()
        metadata = args[0].metadata.data.reset_index()

        facet_args = ['facets_{}'.format(dim) for dim in ['x', 'y']]
        facet_vals = pd.Series({name: kwargs.pop(name, None) for name in facet_args}).dropna()
        
        if facet_vals.size == 0:
            return func(*args, **kwargs)

        if 'groups' in kwargs and kwargs['groups'] not in metagenome.metadata.qual_vars:
            sys.exit('Grouping factor {} not found in metadata.'.format(kwargs['groups']))

        if any(x not in metagenome.metadata.qual_vars for x in facet_vals.tolist()):
            sys.exit('Factor {} not found in metadata.'
                     .format('/'.join(facet_vals.tolist())))

        facet_factors = metadata[facet_vals.tolist()]

        ncol = 1
        if 'facets_y' in facet_vals.index:
            ncol = len(facet_factors[facet_vals['facets_y']].unique())
        
        facet_factors = facet_factors.apply(tuple, axis=1)
        levels = facet_factors.drop_duplicates().sort_values()

        grid = []
        for level in levels:
            samples = metagenome.samples()[facet_factors == level]

            metagenome_subset = metagenome.copy()
            metagenome_subset.subset_samples(sample_names=samples)

            if 'groups' in kwargs:
                metagenome_subset.group_samples(kwargs['groups'])

            plots = func(metagenome_subset, *args[1:], **kwargs)

            facet_suffix = ', '.join(
                ': '.join(info) for info in zip(facet_vals.tolist(), level)
            )

            for p in plots:
                p.title.text = "{} - {}".format(facet_suffix, p.title.text)

            grid += plots

        ncol = max(ncol, len(plots))
        grid = gridplot(grid, ncols=ncol,
                        plot_width=p.plot_width,
                        plot_height=p.plot_height)

        return grid

    return wrapper

def pairwise(func):
    def wrapper(*args, **kwargs):

        (metagenome, factor) = args[:2]
        
        if type(factor) == str:
            factor = metagenome.metadata.data[factor]        

        if 'pairwise' not in kwargs or not kwargs['pairwise']:
            return func(metagenome, factor, *args[2:], **kwargs)

        results = {}
        for (lvl1, lvl2) in combinations(factor.unique(), 2):
            data = metagenome.copy()
            subset = data.samples()[(factor == lvl1) | (factor == lvl2)]
            data.subset_samples(subset)

            name = '{}-vs-{}'.format(lvl1, lvl2)
            results[name] = func(data, factor.loc[subset], *args[2:], **kwargs)

        return pd.DataFrame(results)

    return wrapper
