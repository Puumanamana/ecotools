import numpy as np
import pandas as pd
import h5py

from ecotools.core.biotable import BioTable
from ecotools.util import get_str_dtype


class MetadataTable(BioTable):
    def __init__(self, *qual_vars, **kwargs):

        BioTable.__init__(self, **kwargs)
        self.names = ['samples', 'covariates']
        
        self.quant_vars = self.data.select_dtypes(include='number').columns.to_numpy()
        self.qual_vars = self.data.drop(self.quant_vars, axis=1).columns.to_numpy()
        
    def __repr__(self):
        return "Metadata table: {} samples, {} columns ({} factors, {} covariates)".format(
            *self.data.shape, len(self.qual_vars), len(self.quant_vars)
        )
    
    def to_h5(self, output='metadata.h5'):
        qual_dtype = get_str_dtype(self.data[self.qual_vars])
        col_dtype = get_str_dtype(self.data.columns)
        sample_dtype = get_str_dtype(self.data.index)

        with h5py.File(output, 'a') as handle:
            handle.create_dataset(
                'metadata_qual',
                data=self.data[self.qual_vars].to_numpy().astype(qual_dtype)
            )
            handle.create_dataset(
                'metadata_quant',
                data=self.data[self.quant_vars].to_numpy()
            )            
            handle.create_dataset(
                'qual_vars',
                data=self.qual_vars.astype(col_dtype)
            )
            handle.create_dataset(
                'quant_vars',
                data=self.quant_vars.astype(col_dtype)
            )

            if not 'samples' in handle:
                handle.create_dataset(
                    'samples',
                    data=self.data.index.to_numpy().astype(sample_dtype)
                )

    def group_samples(self, factors):
        # keep the columns with a unique value when aggregated
        qual_vars = np.setdiff1d(self.qual_vars, factors)
        
        cols_to_drop = qual_vars[self.data
                                 .groupby(factors)[qual_vars]
                                 .agg(lambda x: len(set(x)))
                                 .max() > 1]

        self.qual_vars = np.setdiff1d(self.qual_vars, cols_to_drop)

        agg = dict((col, 'mean') if col in self.quant_vars
                   else (col, 'first') for col in self.data.columns)
        
        self.data = self.data.groupby(factors).agg(agg).drop(cols_to_drop, axis=1)
        self.data.index.names = factors

    def factor_data(self, col=None):
        if col is None:
            return self.data[self.qual_vars]
        return self.data[col]

    def subset_rows(self, x):
        self.data = self.data.loc[x]
        self.data = self.data.loc[:, self.data.count() > 0]

    def add_var(self, name, values):
        try:
            values = pd.to_numeric(values)
            self.quant_vars = np.append(self.quant_vars, name)
        except:
            self.qual_vars = np.append(self.qual_vars, name)

        self.data.loc[:, name] = values

    def drop_vars(self, *names):
        names = list(names)
        self.quant_vars = np.setdiff1d(self.quant_vars, names)
        self.qual_vars = np.setdiff1d(self.qual_vars, names)
        self.data.drop(names, axis=1, inplace=True)
