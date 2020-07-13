import numpy as np
from bokeh.models import Arrow, VeeHead, Label


def arrows(pos_i=None, pos_f=None, p=None, names=None, color='black'):

    if pos_i is None:
        pos_i = np.zeros_like(pos_f)
    
    for i, (x0, y0, xf, yf) in enumerate(np.hstack((pos_i, pos_f))):
        p.add_layout(
            Arrow(end=VeeHead(fill_color=color, line_color=color, size=10),
                  x_start=x0, y_start=y0, x_end=xf, y_end=yf)
        )

        if names is not None:
            vec_norm = np.sqrt((xf - x0)**2 + (yf - y0)**2)
            
            labels = Label(x=xf, y=yf,
                           text=names[i],
                           text_color=color,
                           x_offset=vec_norm*0.2, y_offset=vec_norm*0.2)

            p.add_layout(labels)

    return p
