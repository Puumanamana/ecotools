import sys

import numpy as np
import pandas as pd
import h5py

from ecotools.core.biotable import BioTable
from ecotools.util import get_str_dtype


class TaxonomyTable(BioTable):
    ranks = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    
    def __init__(self, species_path='', **kwargs):
        BioTable.__init__(self, **kwargs)
        self.names = ['OTUs', 'levels']

    def _to_h5(self, output):
        
        dtype = get_str_dtype(self.data)
        otu_dtype = get_str_dtype(self.data.index)

        with h5py.File(output, 'a') as handle:
            handle.create_dataset(
                'taxonomy',
                data=self.data.to_numpy().astype(dtype)
            )
            handle.create_dataset(
                'tax_levels',
                data=self.data.columns.to_numpy().astype('S8')
            )

            if not 'OTUs' in handle:
                handle.create_dataset(
                    'OTUs',
                    data=self.data.index.to_numpy().astype(otu_dtype)
                )

    def clean_labels(self, inplace=True, trim=False, maxlen=20):
        data = self.data.replace(to_replace=r'[\-\(\)/]', value='_', regex=True)

        if trim:
            data = (data.fillna('')
                    .applymap(lambda x: '{}{}'.format(x[:maxlen], '...'*(len(x)>maxlen))))

        if not inplace:
            return data

        self.data = data

    def get_ranks(self, info):
        '''
        Get either of the rank, values in info (logical or)
        '''

        info = [(rank.title(), pd.Series(vals).str.lower()) for (rank, vals) in info]
        data_lower = self.data.fillna('').applymap(lambda x: x.lower())
        
        conditions = [
            data_lower[rank].isin(vals)
            for (rank, vals) in info
        ]

        if len(conditions) == 1:
            condition = conditions[0]
        else:
            condition = np.logical_or(*conditions)

        return self.data.index[condition]

    def get_clade(self, other, threshold=0.5):
        if 'Species' not in self.columns:
            sys.exit('Error: No species information')

        species = self.data.Species.dropna().str.split('/').explode()
        species = self.data.reindex(index=species.index).Genus + ' ' + species

        is_clade = (species.str.lower()
                    .isin(other.str.lower())
                    .groupby(level=0).agg('mean'))
        
        clade_otus = is_clade.index[is_clade > threshold]
        
        hits = self.data.Species.loc[clade_otus]

        return hits

    def get_species(self, other):
        return self.get_clade(other, threshold=0)
