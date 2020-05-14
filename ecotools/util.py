import errno
import os
from math import gcd
from pathlib import Path

import numpy as np
from sklearn.mixture import GaussianMixture

from bokeh.palettes import linear_palette, Turbo256

def get_palette(n, random=False):
    palette = linear_palette(Turbo256, n)

    if random:
        p = next(i for i in range(n//10, n) if gcd(i, n) == 1)
        return [palette[i % n] for i in range(0, p*n, p)]

    return palette
 
def get_attributes(obj, keyword):
    return [x for x in dir(obj) if keyword.lower() in x.lower()]

def find_pipeline_files(run_dir, otu_thresh=100):
    files = {
        'abd_path': Path(run_dir, 'Results/main/details', f'abundance_table_{otu_thresh}.shared'),
        'tax_path': Path(run_dir, 'Results/main/details', f'annotations_{otu_thresh}.taxonomy'),
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

    

    
