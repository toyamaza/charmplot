from setuptools import setup, find_packages

setup(name='charmplot',
      version='1.0',
      author='Miha Muskinja',
      author_email='miha.muskinja@cern.ch',
      packages=find_packages(),
      scripts=[
          'charmplot/scripts/compare_multijet.py',
          'charmplot/scripts/get_fake_rates.py',
          'charmplot/scripts/get_production_fractions.py',
          'charmplot/scripts/plot_data_mc.py',
          'charmplot/scripts/plot_mc_mc.py',
          'charmplot/scripts/plot_mc.py',
          'hfroofit/scripts/likelihood_fit.py',
      ],
      install_requires=[
          'PyYAML',
          'lhapdf',
      ]
      )
