import sys

import pandas as pd
import numpy as np

from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import RandomForestClassifier

from bokeh.plotting import figure
from bokeh.io import save, output_file
from bokeh.layouts import gridplot
from bokeh.palettes import inferno, Category10
from bokeh.models import ColumnDataSource, Arrow, NormalHead, LabelSet
from bokeh.transform import jitter

from load_and_convert import load_meta, parse_config
from factors_vs_16S import parse_args

def boxplot_single(info, cmap, col1, col2=None, scatter_data=None):

    info['IQR'] = info.q3 - info.q1

    info['upper'] = np.min([info['max'], info.q3 + 1.5*info.IQR], axis=0)
    info['lower'] = np.max([info['min'], info.q1 - 1.5*info.IQR], axis=0)

    info['color'] = [cmap[f] for f in info[col1]]

    tooltips = []
    if scatter_data is not None:
        tooltips = scatter_data.drop(['y', col1], axis=1).columns
        tooltips = zip(tooltips, '@'+tooltips)
    
    p = figure(title=f"{col2} vs {col1}",
               background_fill_color="#efefef",
               plot_width=300, plot_height=400,
               tooltips=list(tooltips),
               x_range=info[col1].tolist())

    # stems
    p.segment(col1, 'upper', col1, 'q3', line_color="black", source=info)
    p.segment(col1, 'lower', col1, 'q1', line_color="black", source=info)

    # boxes
    p.vbar(x=col1, width=0.7, bottom='q2', top='q3',
           line_color="black", color='color', source=info)
    p.vbar(x=col1, width=0.7, bottom='q1', top='q2',
           line_color="black", color='color', source=info)

    if scatter_data is not None:
        scatter_data = scatter_data.sample(frac=1).groupby(col1).head(500)
        scatter_data['color'] = [cmap[x] for x in scatter_data[col1]]

        p.circle(x=jitter(col1, 0.2, range=p.x_range), y='y',
                 line_color='black', fill_color='color', alpha=0.5,
                 source=scatter_data)

    # # whiskers (almost-0 height rects simpler than segments)
    # h = np.abs(info.q3).mean()
    # p.rect(x=col1, y='lower', width=0.2, height=0.01*h, color="black", source=info)
    # p.rect(x=col1, y='upper', width=0.2, height=0.01*h, color="black", source=info)

    p.xaxis.major_label_orientation = "vertical"

    return p

def bokeh_variation_with_factor(chem, meta, factor_name, fig_dir=None):
    descr = (chem.groupby(meta[factor_name])
             .describe()
             .rename(columns={'25%': 'q1', '50%': 'q2', '75%': 'q3'}))

    if len(descr) > 11:
        colors = inferno(n=len(descr))
    else:
        colors = Category10[len(descr)]

    cmap = {f: colors[i] for (i, f) in enumerate(descr.index)}

    plots = []

    for col in chem.columns:
        data = descr[col].dropna().reset_index()

        if data['std'].max() == 0 or data[['q1','q2','q3']].max().max() == 0:
            continue
        if len(data) > 1:
            scatter_data = meta.assign(y=chem[col]).dropna()
            p = boxplot_single(data, cmap, factor_name, col, scatter_data=scatter_data)
            plots.append(p)

    grid = gridplot(plots, ncols=4)

    output_file(f"{fig_dir}/chem_with_{factor_name}.html")
    save(grid)

def apply_decomp(model_type, X, y=None, fig_dir=None):

    X = (X - X.mean()) / X.std()
    
    if model_type.upper() == 'LDA':
        model = LinearDiscriminantAnalysis(n_components=2)
        model.fit(X, y)
        coeffs =  model.transform(np.identity(X.shape[1]))
    elif model_type.upper() == 'PCA':
        model = PCA(n_components=2)
        model.fit(X)
        coeffs = model.components_.T
    else:
        sys.exit(f"Unknown model: {model_type}")

    components_df = pd.DataFrame(
        model.transform(X),
        columns=[model_type[:2]+'1', model_type[:2]+'2'])
    
    coeffs_df = pd.DataFrame(coeffs, columns=['x', 'y'], index=X.columns)

    biplot(components_df.copy(), coeffs_df, model.explained_variance_ratio_,
           decomp=model_type,
           meta=y,
           fig_dir=fig_dir)
        
def make_classifier(model_type, X, y):
    if model_type.upper() == 'RF':
        model = RandomForestClassifier(n_estimators=100)
        model.fit(X, y)
        coeffs = pd.Series(model.feature_importances_, index=X.columns)
    else:
        sys.exit(f"Unknown model: {model_type}")

    print(coeffs.sort_values(ascending=False))

def biplot(score, coeff, expl_var, decomp='PCA', meta=None, fig_dir=None):

    coeff /= (np.abs(coeff).max() - np.abs(coeff).min()) * 1.5
    score /= (score.max()-score.min())

    coeff['norm'] = coeff.x**2 + coeff.y**2
    coeff = coeff.sort_values(by='norm', ascending=False).iloc[:8]
    
    labels = meta.unique()

    if len(labels) < 11:
        colors = Category10[len(labels)]
    else:
        colors = inferno(len(labels))
    color_map = {label: colors[i] for i, label in enumerate(labels)}

    score[meta.name] = meta.values
    score['color'] = [color_map[label] for label in meta.values]

    p = figure(plot_width=800, plot_height=800, x_range=(-1, 1), y_range=(-1, 1),
               title='Covariate contribution to the first 2 components')
    p.circle(x=score.columns[0], y=score.columns[1], color='color', alpha=0.25, size=5,
             legend_group=meta.name, source=score)
    
    for (x, y, _) in coeff.values:
        p.add_layout(Arrow(end=NormalHead(fill_color="orange", size=12, line_width=1),
                           x_start=0, y_start=0,
                           x_end=x, y_end=y))

    labels = LabelSet(x='x', y='y', text='index', level='glyph',
                      x_offset=5, y_offset=5, render_mode='canvas',
                      source=ColumnDataSource(coeff.reset_index()))

    p.add_layout(labels)
    p.xaxis.axis_label = 'Component 1 ({:.1%})'.format(expl_var[0])
    p.yaxis.axis_label = 'Component 2 ({:.1%})'.format(expl_var[1])

    output_file(f"{fig_dir}/biplot_{decomp}.html")
    save(p)

def bokeh_decomp(chem, factor, fig_dir=None):

    cols = chem.columns[chem.count() > 400]
    idx = chem[cols].dropna(how='any', axis=0).index

    print("Chemistry shape: {} samples, {} variables".format(len(idx), len(cols)))
    
    make_classifier('RF', chem.loc[idx, cols], factor.loc[idx])
    # apply_decomp('LDA', chem.loc[idx, cols], factor.loc[idx], fig_dir=fig_dir)
    apply_decomp('PCA', chem.loc[idx, cols], factor.loc[idx], fig_dir=fig_dir)

if __name__ == '__main__':
    cfg = parse_config()
    args = parse_args()
    fig_dir = cfg.get('misc', 'fig_dir')

    meta = load_meta(cfg)
    meta = meta[~meta[args.factor].isnull()]
    chem = meta.select_dtypes(['number']).drop('Year', axis=1)
    meta.drop(chem.columns, axis=1, inplace=True)

    bokeh_decomp(chem, meta[args.factor], fig_dir)
    bokeh_variation_with_factor(chem, meta, args.factor, fig_dir)
