from pathlib import Path

def find_pipeline_files(run_dir, otu_thresh=100):
    files = {
        'abd_path': Path(run_dir, 'Results/main/details', f'abundance_table_{otu_thresh}.shared'),
        'tax_path': Path(run_dir, 'Results/main/details', f'annotations_{otu_thresh}.taxonomy'),
        'fasta_path': Path(run_dir, 'Results/main/details', f'otu_repr_{otu_thresh}.fasta'),
        'tree_path': Path(run_dir, 'Results/postprocessing/unifrac', f'otu_repr_{otu_thresh}.tre'),
        'species_path': Path(run_dir, 'Results/postprocessing', f'species_{otu_thresh}.csv')
    }
    
    for (name, path) in files.items():
        if not path.is_file():
            print(f'Warning: Could not find {name} file')

            if name not in {'tree_path', 'species_path'}:
                raise FileNotFoundError
            
    return files
