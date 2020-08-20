import numpy as np
import pandas as pd
import h5py

from ecotools.core.biotable import BioTable
from ecotools.util import get_str_dtype, filter_groups


class MetadataTable(BioTable):
    def __init__(self, *qual_vars, **kwargs):

        BioTable.__init__(self, **kwargs)
        self.names = ['samples', 'covariates']
        
        self.quant_vars = self.data.select_dtypes(include='number').columns.to_numpy()
        self.qual_vars = self.data.drop(self.quant_vars, axis=1).columns.to_numpy()

        for var in self.qual_vars:
            if self.data[var].str.strip().str.isdigit().all():
                max_len = self.data[var].apply(len).max()
                self.data[var] = self.data[var].apply(lambda x: '{:0{}}'.format(int(x), max_len))

        self.data[self.qual_vars] = self.data[self.qual_vars].apply(pd.Categorical)
        
    def __repr__(self):
        return "Metadata table: {} samples, {} columns ({} factors, {} covariates)".format(
            *self.data.shape, len(self.qual_vars), len(self.quant_vars)
        )
    
    def _to_h5(self, output):
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

    def group_samples(self, factors, fn='mean'):
        # keep the columns with a unique value when aggregated
        self.data = (self.data.groupby(factors)
                     .pipe(filter_groups, numeric_data=self.quant_vars, fn=fn))
        self.data.index.names = factors
        
    def factor_data(self, cols=None):
        
        if cols is None:
            cols = self.qual_vars
        if isinstance(cols, str):
            cols = [cols]

        data = self.data[cols]

        return data

    def subset_rows(self, x):
        self.data = self.data.loc[x]
        self.data = self.data.loc[:, self.data.count() > 0]

    def add_var(self, name, values, convert=True, order=None):
        if convert:
            try:
                values = pd.to_numeric(values)
                self.quant_vars = np.append(self.quant_vars, name)
            except:
                values = pd.Categorical(values)
                if order is not None:
                    values.reorder_categories(order, inplace=True)
                self.qual_vars = np.append(self.qual_vars, name)
        else:
            if pd.api.types.is_numeric_dtype(pd.Series(values)):
                self.quant_vars = np.append(self.quant_vars, name)
            else:
                self.qual_vars = np.append(self.qual_vars, name)

        self.data.loc[:, name] = values

    def drop_vars(self, *names):
        names = list(names)
        self.quant_vars = np.setdiff1d(self.quant_vars, names)
        self.qual_vars = np.setdiff1d(self.qual_vars, names)
        self.data.drop(names, axis=1, inplace=True)
