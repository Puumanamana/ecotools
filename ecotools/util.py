import errno
import warnings
import os
from pathlib import Path
from math import atan2, pi, gcd

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import scipy.stats
import scipy.optimize
from scipy.spatial import ConvexHull

from bokeh.plotting import figure
from bokeh.io import save, output_file
from bokeh.models import tickers
import bokeh.palettes

def get_str_dtype(data):

    def str_len(x): return len(str(x))
    
    if len(data.shape) == 1:
        max_len = data.map(str_len).max()
    else:
        max_len = data.applymap(str_len).max().max()

    return 'S{}'.format(max_len)

def get_palette(n, palette='Turbo256', random=False, paired=False):

    if paired and n > 20:
        print('Cannot pick a paired palette: Too many colors for `Category20`')

    if paired and n <= 20:
        colors = getattr(bokeh.palettes, 'Category20')[20][:n]
    elif n <= 10:
        colors = getattr(bokeh.palettes, 'Category10')[10][:n]
    else:
        palette = getattr(bokeh.palettes, palette)
        colors = bokeh.palettes.linear_palette(palette, n)

    if random:
        p = next(i for i in range(n//10, n) if gcd(i, n) == 1)
        return [colors[i % n] for i in range(0, p*n, p)]

    return colors


def filter_groups(grouped, numeric=None, fn='mean', approx=False):

    if approx:
        # Just look at the largest group
        is_uniq = grouped.get_group(grouped.size().idxmax()).nunique() == 1
        cols_to_keep = grouped.obj.columns[is_uniq].drop(grouped.grouper.names)
    else:
        is_uniq = grouped.nunique().max() == 1
        cols_to_keep = grouped.obj.columns.drop(grouped.grouper.names)[is_uniq]

    data = grouped[cols_to_keep].nth(0)

    if numeric is not None:
        numeric_data = grouped[numeric].agg(fn)
        data = pd.concat([data.drop(columns=numeric, errors='ignore'), numeric_data], axis=1)
        
    return data
    
def get_attributes(obj, keyword):
    return [x for x in dir(obj) if keyword.lower() in x.lower()]


def find_pipeline_files(run_dir, otu_thresh=100):
    files = {
        'abundance_path': Path(run_dir, 'Results/main/details', f'abundance_table_{otu_thresh}.shared'),
        'taxonomy_path': Path(run_dir, 'Results/main/details', f'annotations_{otu_thresh}.taxonomy'),
        'fasta_path': Path(run_dir, 'Results/main/details', f'otu_repr_{otu_thresh}.filter.fasta'),
        'tree_path': Path(run_dir, 'Results/postprocessing/unifrac', f'FastTree_{otu_thresh}.tre'),
        'species_path': Path(run_dir, 'Results/postprocessing', f'species_{otu_thresh}.csv'),
    }
    
    for (name, path) in files.items():
        if not path.is_file():
            print(f'Warning: Could not find {name} file')

            if name not in {'tree_path', 'species_path', 'fasta_path'}:
                raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), str(path))
            
    return files

def guess_subsampling_level(counts, plot=True, min_value=100):

    if np.std(counts) < 10:
        return np.min(counts)
    
    log_counts = np.log10(counts[counts > min_value])
    (_, sd) = scipy.stats.norm.fit(log_counts)

    (hist, bins) = np.histogram(log_counts.values, bins=30)
    i_mode = np.argmax(hist)
    mode = (bins[i_mode] + bins[i_mode+1]) / 2

    level = 10**(mode - 2*sd)

    if plot:
        bokeh_hist(log_counts, estimator=scipy.stats.norm(mode, sd).pdf)

    if level < 1000:
        other_criteria = counts.quantile(0.1, interpolation='nearest')
        level = max(other_criteria, level)
        

    return int(level)

def filter_metagenome(metagenome, inplace=False, relabund=False,
                      taxa_files=None, taxa=None, clade=False,
                      rank=None, abd_thresh=None):

    if not inplace:
        metagenome = metagenome.copy()

    if relabund:
        # Normalize by sample sum
        metagenome.to_relative_abundance()
    
    if taxa_files is not None or taxa is not None:
        metagenome.subset_otus(taxa_files=taxa_files, clade=clade, taxa=taxa)

    if rank is not None:
        # Only show annotation at the `rank` level
        metagenome.group_taxa(rank)

    if abd_thresh is not None:
        otus = metagenome.get_most_abundant_otus(thresh=abd_thresh)
        metagenome.subset_otus(otus)

    if not inplace:
        return metagenome

def fit_hull(X, n=50):
    # Select points that are in the 95th percentile in each direction
    t = np.arange(n) * 2*pi/n
    polar = np.vstack([np.cos(t), np.sin(t)]).T
    proj = np.sum(X*polar[:, None, :], axis=2)
    
    q3 = np.percentile(proj, 75, axis=1)
    iqr = scipy.stats.iqr(proj, axis=1)
    limits = q3 + 1.5*iqr

    excluded = np.any(proj.T > limits, axis=1)

    if sum(~excluded) < 3:
        excluded = np.zeros(X.shape[0], dtype=bool)

    vertices = ConvexHull(X[~excluded, :]).vertices
    hull_pts = X[~excluded][vertices].T

    return hull_pts
    
def fit_ellipse(X, percentile=90):
    
    center = X.mean(axis=0)
    X_center = X - center

    pcs = PCA().fit(X_center).components_
    angle = atan2(pcs[1, 0], pcs[0, 0])

    u = np.array([np.cos(angle), np.sin(angle)])
    v = np.array([-np.sin(angle), np.cos(angle)])

    # scalar product on each axis
    projs = [np.abs((X_center * u).sum(axis=1)), np.abs((X_center * v).sum(axis=1))]
    (a, b) = np.percentile(projs, percentile, axis=1)

    return (center, a, b, angle)


# To remove or improve
def r_fit_ellipse(X, confidence=0.9, n=100):
    from statsmodels.stats.correlation_tools import cov_nearest
    from scipy.stats import f
    
    cov = cov_nearest(np.cov(X.T), method='clipped', threshold=1e-5)
    center = X.mean(axis=0).reshape(-1, 1)
    
    chol_decomp = np.linalg.cholesky(cov)

    radius = np.sqrt(2 * f.ppf(confidence, 2, X.shape[0]-1))

    angles = np.arange(n) * 2 * pi/n
    unit_circle = np.vstack((np.cos(angles), np.sin(angles))).T
    ellipse = (center + radius * (unit_circle.dot(chol_decomp)).T).T

    return ellipse

def bokeh_hist(log_counts, estimator=None):

    hist, edges = np.histogram(log_counts, bins=max(5, len(log_counts)//10), density=False)
    hist_df = pd.DataFrame({'count': hist,
                            'left': edges[:-1],
                            'right': edges[1:]})
    hist_df['interval'] = ['{:,} - {:,}'.format(int(10**left), int(10**right))
                           for left, right in zip(hist_df['left'], hist_df['right'])]

    x_min = int(min(edges))
    x_max = max(4, 1+int(max(edges)))
    
    p = figure(plot_height=800, plot_width=800,
               x_range=[x_min, x_max], tools='hover,box_zoom',
               tooltips=[('Size range', '@interval'),
                         ('#Samples in interval', str('@count'))],
                title='Sample size distribution',
                x_axis_label='Sample read count',
                y_axis_label='Occurrences')

    p.quad(bottom=0, top='count', left='left', 
           right='right', source=hist_df, fill_color='SteelBlue', 
           line_color='black', fill_alpha=0.7,
           hover_fill_alpha=1.0, hover_fill_color='Tan')

    if estimator is not None:
        x = (hist_df.left + hist_df.right) / 2

        peak = hist_df['count'].max()
        y = estimator(x)
        y *= peak / np.max(y)
        p.line(x, y, line_width=3, color='red')

    ticks = list(range(x_min, x_max))
    minor_ticks = np.log10([i*10**j for i in range(1, 10) for j in ticks])
    
    p.xaxis.ticker = tickers.FixedTicker(ticks=ticks, minor_ticks=minor_ticks)
    p.xaxis.major_label_overrides = {tick: '{:,}'.format(int(10**tick)) for tick in ticks}
    p.yaxis.minor_tick_line_color = None

    p.axis.major_label_text_font_size = '12pt'
    p.axis.axis_label_text_font_size = '14pt'
    p.title.text_font_size = '18pt'

    output_file(os.path.expanduser('~/.trash/sample_sizes.html'))
    save(p)

def make_path(outdir, ext, *args, **kwargs):
    info = [str(x) for x in args if x is not None and x]
    info += [f'{k}-{v}' for (k, v) in kwargs.items() if v is not None and v]

    return Path(outdir, '_'.join(info) + ext)
