from math import sin, cos

import numpy as np
import pandas as pd

from ecotools.plotting.facetgrid import BokehFacetGrid
from ecotools.plotting.scatter import swarmplot, scatter
from ecotools.plotting.barplot import stackplot, barplot
from ecotools.plotting.boxplot import boxplot

def sim(N=1000, seed=42):
    df = pd.DataFrame({
        'x': np.random.choice(list("abcdef"), N),
        'y': np.random.randint(0, 100, N),
        'fact': np.random.choice(['alpha', 'beta', 'gamma'], N),
        'col': np.random.choice(['c1', 'c2'], N)
    })

    return df

def test_swarm_1():
    df = sim()
    g = BokehFacetGrid(data=df, hue='fact', width=800)
    g.map(swarmplot, x='x', y='y', tooltips=['fact', 'x'])

def test_swarm_2():
    df = sim()
    g = BokehFacetGrid(data=df, width=800)
    g.map(swarmplot, x='x', y='y')

def test_boxplot_1():
    df = sim()
    g = BokehFacetGrid(data=df, hue='fact', width=800)
    g.map(boxplot, x='x', y='y', tooltips=['fact', 'x'])

def test_boxplot_2():
    df = sim()
    g = BokehFacetGrid(data=df, width=800)
    g.map(boxplot, x='x', y='y')

def test_box_swarm():
    df = sim()
    g = BokehFacetGrid(data=df, hue='fact', width=800)
    g.map(boxplot, x='x', y='y', tooltips=['fact', 'group'])
    g.map(swarmplot, x='x', y='y')

def test_barplot_1():
    df = sim()
    g = BokehFacetGrid(data=df, hue='fact', width=800)
    g.map(barplot, x='x', y='y', tooltips=['fact', 'group'])

def test_barplot_2():
    df = sim()
    g = BokehFacetGrid(data=df, width=800)
    g.map(barplot, x='x', y='y')

def test_stackplot_1():
    df = sim()
    g = BokehFacetGrid(data=df, width=800, hue='fact')
    g.map(stackplot, x='x', y='y')
    
def test_stackplot_2():
    df = sim()
    g = BokehFacetGrid(data=df, width=800, hue='fact')
    g.map(stackplot, x=['x', 'col'], y='y')

def test_scatter_ellipse():
    N = 500
    df = pd.DataFrame(np.random.multivariate_normal([3, 5], [[10,5], [5,4]], N))
    df.columns=['x', 'y']
    df['fact'] = np.random.choice(['a', 'b', 'c'], N, replace=True)
    angle = 0.5
    
    R = np.array([[cos(angle), -sin(angle)], [sin(angle), cos(angle)]])

    df[['x', 'y']] = df[['x', 'y']].dot(R).to_numpy()

    df.loc[df.fact=='a', ['x', 'y']] = df.loc[df.fact=='a', ['x', 'y']].dot(R).to_numpy()
    df.loc[df.fact=='b', ['x', 'y']] = df.loc[df.fact=='b', ['x', 'y']].dot(R).dot(R).to_numpy()
    
    g = BokehFacetGrid(data=df, scale=2, hue='fact')
    g.map(scatter, x='x', y='y', ellipse=True)
    g.save('ellipse_test.html')
    
if __name__ == '__main__':
    test_scatter_ellipse()    
