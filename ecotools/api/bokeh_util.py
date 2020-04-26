import sys

from bokeh.palettes import inferno
from bokeh.models import Legend
from bokeh.io import output_file, save
from bokeh.layouts import gridplot

TOOLS = ['hover', 'box_zoom', 'pan', 'reset']
PADDING = 50
LEG_TXT_LEN = 30
LEG_NROWS = 25
LEG_MAXCOL = 3

def get_palette(n):
    return inferno(n)
 
def get_attributes(obj, keyword):
    return [x for x in dir(obj) if keyword.lower() in x.lower()]

def bokeh_legend_out(func):
    def wrapper(*args, **kwargs):
        p = func(*args, **kwargs)
        
        leg_items = []
        for item in p.legend.items[::-1]:
            label = item.label['value']
            item.label['value'] = label[:LEG_TXT_LEN] + '..'*(len(label) > LEG_TXT_LEN)
            leg_items.append(item)        

        p.legend.items.clear()
        
        col = 0 
        while leg_items:
            if col == LEG_MAXCOL:
                items = [leg_items.pop() for _ in range(len(leg_items))]
            else:
                items = [leg_items.pop() for _ in range(LEG_NROWS) if leg_items]
            legend = Legend(items=items, location=(25, 0))
            p.add_layout(legend, 'right')
            col += 1

        p.legend.label_text_font_size = '8pt'
        p.legend.padding = 0
        p.legend.spacing = 0
        p.legend.border_line_color = None
        # p.legend.label_width = LEG_TXT_LEN

        return p
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
            
    
def bokeh_facets(func):
    '''
    Decorator to make facets with factor levels
    '''

    def wrapper(*args, **kwargs):

        if not type(args[0]).__name__ == 'MetagenomicDS':
            sys.exit('First argument must be of type MetagenomicsDS. Found {}'
                     .format(type(args[0]).__name__))

        metagenome = args[0]

        if 'facets' not in kwargs or kwargs['facets'] is None:
            return func(*args, **kwargs)

        if any(x not in metagenome.metadata.qual_vars for x in kwargs['facets']):
            sys.exit('Facet factor {} not found in metadata.'.format(kwargs['facets']))
            return func(*args, **kwargs)

        facet_factors = list(kwargs.pop('facets'))
        factors = metagenome.metadata.data[facet_factors]
        ncol = 1

        if len(facet_factors) > 1:
            ncol = len(factors[facet_factors[1]].unique())
            factors = factors.apply(tuple, axis=1)
            levels = factors.unique()
        else:
            factors = factors.iloc[:, 0]
            levels = factors.drop_duplicates()

        grid = []
        for level in levels:
            samples = metagenome.samples()[factors == level]

            metagenome_subset = metagenome.copy()
            metagenome_subset.subset_samples(sample_names=samples)
            
            p = func(metagenome_subset, *args[1:], **kwargs)

            if len(facet_factors) == 1:
                facet_suffix = '{}: {}'.format(facet_factors[0], level)
            else:
                facet_suffix = ', '.join(': '.join(info) for info in zip(facet_factors, level))
            p.title.text = "{} - {}".format(facet_suffix, p.title.text)

            grid.append(p)
            
        grid = gridplot(grid, ncols=ncol,
                        plot_width=p.plot_width,
                        plot_height=p.plot_height)

        return grid

    return wrapper

