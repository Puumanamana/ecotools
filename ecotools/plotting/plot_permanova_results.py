from pathlib import Path
import re

import numpy as np
import pandas as pd

from bokeh.io import output_file, save
from bokeh.plotting import figure
from bokeh.layouts import gridplot
from bokeh.models import FactorRange, Band, ColumnDataSource
from bokeh.transform import factor_cmap
from bokeh.palettes import Dark2, Accent, Category10

from ecotools.load_and_convert import parse_config
from ecotools.factors_vs_16S import parse_args

def load_results(folder, factor):
    data = []

    for f in folder.glob(f"adonis_{factor}_*.csv"):
        df = pd.read_csv(f, index_col=0)
        (lvl1, lvl2) = re.findall("_([^_]*\-vs\-[^_]*)$", f.stem)[0].split('-vs-')
        df['comparison'] = "{}-vs-{}".format(lvl1, lvl2)
        df['level_1'] = lvl1
        df['level_2'] = lvl2
        df['neg_log_padj'] = -np.log10(df['p.adjusted'])

        data.append(df)

    data = pd.concat(data)
    data.index.name = 'season'
    data = data.reset_index().set_index(['season', 'comparison'])

    full_index = pd.MultiIndex.from_product(data.index.levels)
    
    return data.reindex(full_index)


