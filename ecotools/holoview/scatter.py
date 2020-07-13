import holoviews as hv
from holoviews import opts
hv.extension('bokeh')

def scatter(data, x, y, hue=None):
    key_dimensions   = [(xi, xi) for xi in [x, hue]]
    value_dimensions = [(y, y)]

    macro = hv.Table(data, key_dimensions, value_dimensions)

    scatter = macro.to.scatter(x, y).overlay(hue)
    
    scatter.opts(
        opts.Scatter(color=hv.Cycle('Category20'), line_color='k', size=10,
                     show_grid=True, width=700, height=400),
        opts.NdOverlay(legend_position='left', show_frame=False))

    return scatter

    
