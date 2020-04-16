from setuptools import setup, find_packages

setup(name='charmplot',
      version='1.0',
      author='Miha Muskinja',
      author_email='miha.muskinja@cern.ch',
      packages=find_packages(),
      scripts=[
          'charmplot/scripts/charm_plot.py',
          'charmplot/scripts/plot_mc_mc.py',
          'charmplot/scripts/plot_mc.py',
      ],
      )
