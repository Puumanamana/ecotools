import argparse

import pandas as pd
from scipy.stats import pearsonr

import matplotlib.pyplot as plt

from ecotools.load_and_convert import load_h5, parse_config
from ecotools.util import group_by_rank

TOOLS = ['hover', 'box_zoom', 'reset']
OUTDIR = '../outputs'

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--chem', type=str, default='Si')
    args = parser.parse_args()

    return args

def get_correlated_otus(shared, tax, chem, top=4):

    X = {}
    corrs = []

    for lvl in tax.columns:
        print("Processing", lvl)
        X[lvl] = group_by_rank(shared, tax, lvl)[0]

        for i, otu in enumerate(X[lvl].index):
            otu_abd = X[lvl].loc[otu]

            nz = otu_abd > 0
            if sum(nz) < 20:
                continue
            score = pearsonr(otu_abd, chem)[0]
            corrs.append([otu, score, lvl, abs(score)])

    corrs = pd.DataFrame(corrs, columns=['name', 'score', 'level', 'abs_score'])
    corrs = (corrs.groupby('level')
             .apply(lambda x: x.nlargest(top, 'abs_score'))
             .drop(['abs_score', 'level'], axis=1))
    print(corrs)

    fig, ax = plt.subplots(top, tax.shape[1], figsize=(20, 15))
    for i, lvl in enumerate(tax.columns):
        for j, (name, score) in enumerate(corrs.loc[lvl].to_numpy()):
            ax[j, i].scatter(chem, X[lvl].loc[name], s=10, alpha=0.7)
            ax[j, i].set_title('Correlation = {:.2f}'.format(score))
            ax[j, i].set_ylabel(name)

    plt.suptitle("Top {} OTUs correlated with {} (0 removed)".format(top, chem.name))
    plt.subplots_adjust(hspace=0.5, wspace=0.5)
    plt.show()

if __name__ == '__main__':
    args = parse_args()
    cfg = parse_config()

    thresh = cfg.getint('misc', 'otu_thresh')
    
    (shared, taxonomy, _, metadata) = load_h5(cfg.get('misc', 'h5'), norm_fn="assr")

    # Most abundant chemical data: Si, Calcium, Potassium, Sodium, Sulfate, Magnesium, Fluoride, Chloride, Sr, V
    print(metadata.columns)

    meta = metadata[args.chem].dropna().astype(float)
    
    shared_subset = get_correlated_otus(shared.loc[meta.index], taxonomy.drop('Kingdom', axis=1), meta)

    
