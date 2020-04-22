import numpy as np
from bokeh.palettes import linear_palette, Turbo256
from bokeh.models import Legend
from bokeh.io import output_file, save

TOOLS = ['hover', 'box_zoom', 'reset']

def get_palette(n):    
    palette = linear_palette(Turbo256, n)
    palette = np.random.choice(palette, n, replace=False).tolist()

    return palette

def bokeh_legend_out(func):
    def wrapper(*args, **kwargs):
        p = func(*args, **kwargs)
        
        leg_items = []
        for item in p.legend.items:
            label = item.label['value']

            if len(label) > 50:
                item.label['value'] = label[:50] + '...'

            leg_items.append(item)            


        # Not quite right, need to work on the formula
        glyph_height = 20
        offset = 0
        # offset = (p.plot_height - len(leg_items)*p.legend.label_height) / 2

        legend = Legend(items=leg_items, glyph_height=glyph_height, location=(20, -offset))
        p.legend.items.clear()
        p.add_layout(legend, 'right')
        p.legend.label_text_font_size = '8pt'
        
        return p
    return wrapper

def bokeh_save(func):
    def wrapper(*args, **kwargs):
        p = func(*args, **kwargs)

        if 'output' in kwargs and kwargs['output'] is not None:
            output_file(kwargs['output'])
            save(p)
        else:
            return p

    return wrapper
            
        
    
