import sys
import pandas as pd
import h5py

import sklearn.preprocessing
from skbio import diversity
from skbio.stats import subsample_counts

from ecotools.core.biotable import BioTable

def get_str_dtype(data):

    def str_len(x): return len(str(x))
    
    if len(data.shape) == 1:
        max_len = data.map(str_len).max()
    else:
        max_len = data.applymap(str_len).max().max()

    return 'S{}'.format(max_len)


class AbundanceTable(BioTable):

    def __init__(self, **kwargs):
        BioTable.__init__(self, **kwargs)
        
        self.names = ['samples', 'OTUs']
        self.alpha_diversity = None
        self.distance_matrix = None
        self.raw_sample_sizes = self.data.sum(axis=1)        

    def _to_h5(self, output):

        sample_dtype = get_str_dtype(self.data.index)
        otu_dtype = get_str_dtype(self.data.columns)

        with h5py.File(output, 'a') as handle:
            handle.create_dataset(
                'abundance',
                data=self.data.to_numpy()
            )

            if not 'samples' in handle:
                handle.create_dataset(
                    'samples',
                    data=self.data.index.to_numpy().astype(sample_dtype)
                )

            if not 'OTUs' in handle:
                handle.create_dataset(
                    'OTUs',
                    data=self.data.columns.to_numpy().astype(otu_dtype)
                )

    def sample_sizes(self):
        return self.data.sum(axis=1)

    def subset_rows(self, x):
        self.data = self.data.loc[x]
        self.data = self.data.loc[:, self.data.sum() > 0]
        
        try:
            self.raw_sample_sizes = self.raw_sample_sizes.loc[self.data.index]
        except KeyError:
            self.raw_sample_sizes = None
    
    def subset_cols(self, x):
        self.data = self.data.loc[:, x]
        self.data = self.data.loc[self.data.sum(axis=1) > 0, :]

    def to_relative_abundance(self):
        same_idx = self.data.index.isin(self.raw_sample_sizes.index)
        
        if sum(same_idx) == 0:
            print('Cannot normalize by sample size: data already aggregated')
            return
        self.data = (self.data.T / self.raw_sample_sizes.loc[self.index]).T
        
    def normalize(self, method, axis='features'):
        data = self.data
        if axis == 'samples':
            data = data.T

        if not method:
            return
        elif method == 'wisconsin':
            data = data / data.max(axis=0)
            data = (data.T / data.sum(axis=1)).T
        else:
            normalizer = getattr(sklearn.preprocessing, method)()
            data = pd.DataFrame(normalizer.fit_transform(data),
                                index=data.index,
                                columns=data.columns)

        if axis == 'samples':
            data = data.T

        self.data = data

    def subsample(self, level):
        dropped = []

        for (i, row) in enumerate(self.data.to_numpy()):
            try:
                row_subsampled = subsample_counts(row, level)
            except ValueError:
                dropped.append(i)
                continue
                
            self.data.iloc[i] = row_subsampled

        self.data.drop(self.data.index[dropped], inplace=True)

    def get_most_abundant_otus(self, thresh=0.01):
        proportions = (self.data.T / self.data.sum(axis=1)).max(axis=1)
        main_otus = proportions.sort_values(ascending=False).index[proportions > thresh]

        return main_otus

    def get_most_variable_otus(self, n=100):
        proportions = self.data / self.data.sum()
        variance = proportions.std().sort_values(ascending=False)

        return variance.index[:n]

    def select_otus(self, criteria='prevalence', n=100):
        if criteria == 'prevalence':
            scores = (self.data > 0).sum()
        elif criteria == 'abundance':
            scores = self.data.sum()
        elif criteria == 'variance':
            scores = self.data.std()
        else:
            sys.exit('Unknown criteria')

        scores.sort_values(ascending=False, inplace=True)
        return scores.index[:n]


    def group_samples(self, groups, fn='mean'):
        self.data = self.data.groupby(groups).agg(fn)

        if str(fn) == 'sum' or (callable(fn) and fn.__name__ == 'sum'):
            self.raw_sample_sizes = self.raw_sample_sizes.groupby(groups).agg(fn)
    
    def compute_alpha_diversity(self, metric):

        if metric == 'total_abundance':
            self.alpha_diversity = self.data.sum(axis=1)
        elif metric == 'richness':
            self.alpha_diversity = (self.data > 0).sum(axis=1)
        else:
            self.alpha_diversity = self.data.apply(getattr(diversity.alpha, metric), axis=1)

        self.alpha_diversity = pd.Series(self.alpha_diversity, index=self.index)
        self.alpha_diversity.name = metric
