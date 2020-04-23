import numpy as np
from bokeh.palettes import linear_palette, Turbo256
from bokeh.models import Legend
from bokeh.io import output_file, save

TOOLS = ['hover', 'box_zoom', 'reset']
PADDING = 200
LEG_TXT_LEN = 20
LEG_NROWS = 20
LEG_MAXCOL = 3

def get_palette(n):    
    palette = linear_palette(Turbo256, n)
    palette = np.random.choice(palette, n, replace=False).tolist()

    return palette

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
        # p.legend.label_text_line_height = 0.1
        p.legend.padding = 0
        p.legend.spacing = 0
        # p.legend.glyph_height = 5
        p.legend.border_line_color = None
        # p.legend.label_ = 'center'
        # p.legend.label_width = LEG_TXT_LEN

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
            
        
    
