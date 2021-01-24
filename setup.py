from setuptools import setup, find_packages

setup(name='charmplot',
      version='1.0',
      author='Miha Muskinja',
      author_email='miha.muskinja@cern.ch',
      packages=find_packages(),
      scripts=[
          'charmplot/scripts/compare_multijet.py',
          'charmplot/scripts/extract_error_bands_from_rivet_dat_files.py',
          'charmplot/scripts/extract_production_fractions_from_rivet_dat_files.py',
          'charmplot/scripts/extract_variations_from_rivet_dat_files.py',
          'charmplot/scripts/get_fake_rates.py',
          'charmplot/scripts/get_post_fit_yields.py',
          'charmplot/scripts/get_production_fractions.py',
          'charmplot/scripts/get_simple_xsec.py',
          'charmplot/scripts/get_wplusd_xsec.py',
          'charmplot/scripts/plot_data_mc.py',
          'charmplot/scripts/plot_mc_mc.py',
          'charmplot/scripts/plot_mc.py',
          'charmplot/scripts/plot_post_fit.py',
          'gen/generate_wplusd_config.py',
          'hfroofit/scripts/likelihood_fit.py',
      ],
      install_requires=[
          'PyYAML',
      ]
      )
