import os
import pandas as pd
import numpy as np

DB_FILE = os.path.expanduser("~/databases/pathogens/PATRIC_website/human-pathogens_lineage.csv")

def get_pathogens(sp_file):
    species = pd.read_csv(sp_file, index_col=0).Species.dropna()
    db = pd.read_csv(DB_FILE)['species']

    patho_otus = set()
    for i, sp in enumerate(species.to_numpy()):
        res = np.mean([db.str.contains(hit, case=False).sum() > 0 for hit in sp.split('/')])

        if res > 0.5:
            patho_otus.add(i)
            print(len(patho_otus), end='\r')

    print("{} / {} potential pathogenic OTUs".format(len(patho_otus), len(species)))
    return species.iloc[list(patho_otus)].sort_index()
