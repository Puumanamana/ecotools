import sys

sys.path.append('..')
from base_data_classes import Data, MetadataTable, AbundanceTable, TaxonomyTable
from metagenomics_ds import MetagenomicDS

def test_data():
    data = Data(ds='test', path='abc')

    assert data.ds == 'test'
    assert data.shape == None
    assert not data.valid_input()

def test_metadata():
    metadata = MetadataTable(path='test_meta.csv', outdir='/tmp')
    metadata.to_h5()
    metadata = MetadataTable(outdir='/tmp')

    metadata.h5.unlink()

    assert metadata.data.shape == (10, 4)

def test_abd_data():
    abd_data = AbundanceTable(path='test.shared', outdir='/tmp')
    abd_data.to_h5()
    abd_data = AbundanceTable(outdir='/tmp')

    abd_data.h5.unlink()

    assert abd_data.data.shape == (11, 18)

def test_tax_data():
    tax_data = TaxonomyTable(path='test.taxonomy', outdir='/tmp')
    tax_data.to_h5()
    tax_data = TaxonomyTable(outdir='/tmp')

    tax_data.h5.unlink()
    assert tax_data.data.shape == (17, 6)

def test_metagenomic_data():
    metagenome = MetagenomicDS(abd_path='test.shared',
                               tax_path='test.taxonomy',
                               meta_path='test_meta.csv',
                               outdir='/tmp')
    metagenome.to_h5()
    metagenome = MetagenomicDS(outdir='/tmp')
    metagenome.h5.unlink()

    shapes = metagenome.shapes()

    assert shapes['otus'] == 17
    assert shapes['samples'] == 10
    assert shapes['metadata'] == 4

if __name__ == '__main__':
    test_metagenomic_data()
