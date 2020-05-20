from pathlib import Path

class BioData:
    def __init__(self, name='data', data=None, output=None):
        self.data = data
        self.outdir = './'

    def __getattr__(self, key):
        if key in dir(self):
            return getattr(self, key)

        if not key.startswith('__') and 'data' in dir(self):
            data_obj = getattr(self, 'data')

            if key in dir(data_obj):
                return getattr(data_obj, key)

        raise AttributeError(key)

class BioTable(BioData):
    def __init__(self, **kwargs):
        BioData.__init__(self, **kwargs)
        self.names = ['samples', 'OTUs']
        self.title = ''

    def __repr__(self):
        info = [' '.join(map(str, x)) for x in zip(self.data.shape, self.names)]
        return ('{}: {} x {}'.format(self.__class__.__name__, *info))

    def get_column_format(self, index_name, col_name, value_name):
        data = self.data.stack().rename(value_name)
        data.index.names = [index_name, col_name]
        return data

    def subset_rows(self, x):
        self.data = self.data.loc[x]

    def subset_cols(self, x):
        self.data = self.data.loc[:, x]        

    def to_csv(self, path):
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        self.data.to_csv(path)

