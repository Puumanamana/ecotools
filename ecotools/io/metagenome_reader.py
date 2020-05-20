from pathlib import Path

from ecotools.core.metagenomics_ds import MetagenomicDS
from ecotools.io.txt import from_txt
from ecotools.io.binary import from_h5


def metagenome_from_files(
        outdir='.', ds='project',
        abundance_path=None,
        taxonomy_path=None,
        metadata_path=None,
        fasta_path=None,
        **kwargs):

    inputs = {}
    
    for (name, path) in {'abundance': abundance_path,
                         'taxonomy': taxonomy_path,
                         'metadata': metadata_path,
                         'fasta': fasta_path}.items():
        
        if path is None:
            continue

        if Path(path).suffix == '.h5':
            data = from_h5(path, kind=name)
        else:
            data = from_txt(path, **kwargs)

        inputs[name] = data

    return MetagenomicDS(ds=ds, outdir=outdir, **inputs)
