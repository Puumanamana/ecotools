from setuptools import setup, find_packages

setup(name='ecotools',
      version='0.0',
      description='Tool for ecological data visualization and analysis',
      long_description=open('README.rst').read(),
      url='https://github.com/Puumanamana/ecotools',
      packages=find_packages(),
      entry_points={'console_scripts': ['ecotools=ecotools.main:main']},
      license='Apache License 2.0',
      zip_safe=False,
      python_requires='>=3.6',
      test_requires=['pytest', 'pytest-cov'],
      install_requires=[
          'numpy',
          'pandas',
          'h5py',
          'argparse',
          'bokeh>=1.4',
          'scipy',
          'scikit-learn',
          'scikit-bio',
          'rpy2>=3.3.0',
          'tzlocal',
          'pyyaml',
      ])
