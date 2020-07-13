from pathlib import Path

from ecotools.core.metagenomics_ds import MetagenomicDS
from ecotools.io.txt import from_txt
from ecotools.io.binary import from_h5


def metagenome_from_files(
        outdir='.', h5dir='.',
        ds='project',
        abundance_path=None,
        taxonomy_path=None,
        species_path=None,
        metadata_path=None,
        fasta_path=None,
        tree_path=None,
        parser_args={},
        **kwargs):

    inputs = {}

    for (name, path) in {'abundance': abundance_path,
                         'taxonomy': taxonomy_path,
                         'metadata': metadata_path,
                         'sequences': fasta_path,
                         'tree': tree_path}.items():

        h5_path = Path(outdir, f'{ds}_h5-data', f'{name}.h5')

        if h5_path.is_file():
             data = from_h5(h5_path, kind=name)
        elif path is None:
            continue
        else:
            data = from_txt(path, species_path=species_path, **parser_args)

        inputs[name] = data

    mg = MetagenomicDS(ds=ds, outdir=outdir, **inputs, **kwargs)

    return mg
