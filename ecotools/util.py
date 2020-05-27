import errno
import os
from math import gcd
from pathlib import Path
from math import atan2

import numpy as np
from sklearn.mixture import GaussianMixture
from sklearn.decomposition import PCA

import bokeh.palettes

def get_str_dtype(data):

    def str_len(x): return len(str(x))
    
    if len(data.shape) == 1:
        max_len = data.map(str_len).max()
    else:
        max_len = data.applymap(str_len).max().max()

    return 'S{}'.format(max_len)

def get_palette(n, palette='Turbo256', random=False):
    palette = getattr(bokeh.palettes, palette)
    colors = bokeh.palettes.linear_palette(palette, n)

    if random:
        p = next(i for i in range(n//10, n) if gcd(i, n) == 1)
        return [colors[i % n] for i in range(0, p*n, p)]

    return colors

def elt_or_nothing(l):
    if len(set(l)) == 1:
        return l.iloc[0]
    return None
 
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

def guess_subsampling_level(sample_sizes):
    model = GaussianMixture(n_components=2)
    model.fit(np.log10(sample_sizes.to_numpy().reshape(-1, 1)))

    i_max = np.argmax(model.means_)
    level = 10**(model.means_[i_max] - 3*model.covariances_[i_max, 0, 0])

    if level < 2000:
        level = max(sample_sizes.quantile(0.1), level)

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

def fit_ellipse(X):
    
    mu = X.mean(axis=0)
    X_center = X - mu

    pcs = PCA().fit(X_center).components_
    angle = atan2(pcs[1, 0], pcs[0, 0])

    u = np.array([np.cos(angle), np.sin(angle)])
    v = np.array([-np.sin(angle), np.cos(angle)])

    # scalar product on each axis
    projs = [np.abs((X_center * u).sum(axis=1)), np.abs((X_center * v).sum(axis=1))]
    (a, b) = np.percentile(projs, 90, axis=1)

    return (mu, a, b, angle)
