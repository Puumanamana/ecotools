from pathlib import Path

import ecotools as ect

SRC_DIR = Path(__file__).resolve().parent
ABD_PATH = Path(SRC_DIR, 'data', 'abundance.shared')
TAX_PATH = Path(SRC_DIR, 'data', 'annotations.taxonomy')
META_PATH = Path(SRC_DIR, 'data', 'metadata.csv')
SP_PATH = Path(SRC_DIR, 'data', 'species.csv')

def test_load_tax():
    data = ect.from_txt(TAX_PATH, species_path=SP_PATH)

    print(data)
    return True

def test_load_shared():
    data = ect.from_txt(ABD_PATH)

    print(data)
    return True

def test_load_meta():
    data = ect.from_txt(META_PATH)

    print(data)
    return True

def test_load_metagenome():
    data = ect.metagenome_from_files(
        abundance_path=ABD_PATH,
        taxonomy_path=TAX_PATH,
        metadata_path=META_PATH,
        species_path=SP_PATH
    )

    print(data)
    return True

def test_from_to_h5():
    data = ect.metagenome_from_files(
        abundance_path=ABD_PATH,
        taxonomy_path=TAX_PATH,
        metadata_path=META_PATH,
        species_path=SP_PATH
    )

    print(data)
    
    data.to_h5()

    data = ect.metagenome_from_files(
        abundance_path='abundance.h5',
        taxonomy_path='taxonomy.h5',
        metadata_path='metadata.h5',
    )

    print(data)
    
    for filename in ['abundance.h5', 'metadata.h5', 'taxonomy.h5']:
        Path(filename).unlink()

    return True

if __name__ == '__main__':
    test_from_to_h5()

