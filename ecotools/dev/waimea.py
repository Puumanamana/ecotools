from pathlib import Path

from metagenomics_ds import MetagenomicDS
from hue_boxplot import diversity_with_meta
from barplot import taxa_barplot

meta_file = '/home/cedric/data/Kiana/Waimea/metadata.csv'
root = '/home/cedric/data/Kiana/Waimea/16S-pipeline_outputs_sub0/Results/main/details'
data = MetagenomicDS(abd_path=f'{root}/abundance_table_100.shared',
                     tax_path=f'{root}/annotations_100.taxonomy',
                     species_path=f'{root}/../../postprocessing/species_100.csv',
                     meta_path=meta_file,
                     qual_vars=['Site'],
                     outdir='/tmp/cedric/ecotools')

print(data)

pathogens = '~/databases/pathogens/PATRIC_website/human-pathogens_lineage.csv'

diversity_with_meta(data, ['Season', 'Site'], norm=True, metric='n_otus', taxa_file=pathogens,
                    output=Path('pathogenic_by-Site_abd-norm-by-sampleSize.html'))
taxa_barplot(data, factor='Site', norm=True, rank='Species', top=20, taxa_file=pathogens,
             output=Path('pathogenic-species_barplot_by-Site_norm-by-sampleSize_top-20.html'))

