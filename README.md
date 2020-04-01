- [charmplot](#charmplot)
- [Installation on Cori](#installation-on-cori)
  - [Install conda](#install-conda)
  - [Create new pyton3 conda environemnt](#create-new-pyton3-conda-environemnt)
  - [Install ROOT](#install-root)
  - [Install charmplot](#install-charmplot)
- [Base analysis script](#base-analysis-script)
- [Examples](#examples)
  - [Basic OS/SS and OS-SS plots](#basic-osss-and-os-ss-plots)
  - [Mass fit](#mass-fit)
- [Stage-out plots to a webpage](#stage-out-plots-to-a-webpage)

## charmplot

Framework for various data analysis and histogram plotting using
[charmpp](https://gitlab.cern.ch/lbnl/CharmPhysics/charmpp) ntuples.

## Installation on Cori

Both python3 and ROOT are required. On Cori, the easiest way of achieving this
is using conda. To install charmplot in conda follow the instructions below.

### Install conda

If you already have a python3 conda environment on Cori you can skip this step
and just activate that environment.

```
module load python
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh .
chmod +x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh
```

### Create new pyton3 conda environemnt

```
conda create --name myenv python=3.7
source activate myenv
```

On all subsequent logins, the above environment can be activated with the following commands:

```
module load python
source activate myenv
```

### Install ROOT

This one may take a while.

```
conda install -c conda-forge root cmake
```

### Install charmplot

Use pip with option `-e` to install the charmplot package.

```
pip install -e <path-to-charmplot>
```

Options `-e` triggers the 'develop' mode where changing the files will have an
immediate effect without re-installing the package.

## Base analysis script

Most plotting and fitting is preformed through the main [charm_plot.py](https://gitlab.cern.ch/lbnl/CharmPhysics/charmplot/-/blob/master/charmplot/scripts/charm_plot.py)
script. The analysis to run is defined with the `-a` command line argument. Input string should match one of the configurations in the
[config](https://gitlab.cern.ch/lbnl/CharmPhysics/charmplot/-/tree/master/charmplot%2Fconfig) folder. Some examples are given below.

## Examples

First copy the example ntuples to a folder in your home directory:

```
mkdir example
cd example
cp /global/cfs/cdirs/atlas/wcharm/example/v4/* .
```

### Basic OS/SS and OS-SS plots

The basic configuration file to run OS/SS and OS-SS plots is [wplusd_powheg.yaml](https://gitlab.cern.ch/lbnl/CharmPhysics/charmplot/-/blob/master/charmplot/config/wplusd_powheg.yaml).
By inspecting the file we see it has 6 channels defined: 2018_el_wplusd_OS, 2018_el_wplusd_SS, 2018_el_wplusd_OS-SS, 2018_mu_wplusd_OS, 2018_mu_wplusd_SS, 2018_mu_wplusd_OS-SS.
Plotted variables are controlled through the linked [charmed_wjets.yaml](https://gitlab.cern.ch/lbnl/CharmPhysics/charmplot/-/blob/master/charmplot/config/variables/charmed_wjets.yaml) file.
Command line argument `-v` instructs the framework to plot only a subset of variables (comma separated). Similarly, `-c` specifies a subset of channels.

To run over example ntuples run the following commands:

```
charm_plot.py -a wplusd_powheg
```

Inspect the created 'wplusd_powheg' folder. It should have 6 subfolders-- one for each channel, containing all variables.

### Mass fit

Configuration file for the mass fit is [wplusd_mass_fit.yaml](https://gitlab.cern.ch/lbnl/CharmPhysics/charmplot/-/blob/master/charmplot/config/wplusd_mass_fit.yaml).
By default, it will perform mass fits both for 'Data - Bkg' and 'Signal MC' in OS/SS and OS-SS regions. Background function choice and parameter limits can
be configured through the configuration file.

To run over example ntuples run the following commands:

```
charm_plot.py -a wplusd_mass_fit -v Dmeson_m
```

<!-- ### QCD Template fit

An example configuration file for the QCD Template fit is [wplusd_fit/electron_pt.yaml](https://gitlab.cern.ch/lbnl/CharmPhysics/charmplot/-/blob/master/charmplot/config/wplusd_fit/electron_pt.yaml).
In this case the 'PT Template' fit is performed for the electron channel. Multiple channels are defined in the configuration file:

- 2018_el_wplusd_OS: OS signal region (standard mT, MET, pT, d0 cuts)
- 2018_el_wplusd_PT_Template_OS: region used to extract the template (no mT cut, inverted MET and d0, relaxed pT)
- 2018_el_wplusd_PT_Loose_OS: region where the likelihood fit is performed (no mT cut, relaxed pT)
- 2018_el_wplusd_PT_Fit_OS: same as above, but used to control the likelihood fit
- 2018_el_wplusd_PT_Extrapolated_OS: OS signal region with extrapolated QCD Multijet background
- 2018_el_wplusd_PTCut_Template_OS: used only to construct QCD Multijet histograms with a pT cut (used for extrapolation in SR)

To run over example ntuples run the following commands:

```
charm_plot.py -a wplusd_fit/electron_pt -v "Dmeson_m,lep_pt,met_met"
``` -->

## Stage-out plots to a webpage

charm_plot.py contains an option to automatically stage-out plots to a webpage.
Run charm_plot.py with the `--stage-out` command line argument to trigger this.

https://portal.nersc.gov/project/atlas/wcharm/

