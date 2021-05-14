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

## Example D+ Plots

Copy the charmpp output in a clean folder:
```
mkdir <some_dit>
cd <some_dir>
cp /global/cfs/cdirs/atlas/wcharm/charmpp_output/Dplus_2021_05_14/*.root .
```

_NOTE:_ the `--split-signal-samples` option is needed in the examples because the samples were processed with the `bin_in_truth_D_pt: true` option in charmpp. This splits the `Matched` category into `Matched_truth_pt_binX` and the inclusive one is not explicitly saved.

### Basic data / MC plots with truth matching

```
generate_wplusd_config.py -a wplusd -s truth --extra-rebin 2.5 --split-signal-samples
plot_data_mc.py -a wplusd.yaml -v Dmeson_m --nology
```

### Basic data / MC plots without truth matching

```
generate_wplusd_config.py -a wplusd_flavor -s flavor --extra-rebin 2.5 --split-signal-samples
plot_data_mc.py -a wplusd_flavor.yaml -v Dmeson_m --nology
```

### SPG comparison

```
generate_wplusd_config.py -a dplus_spg -s spg_comparison --extra-rebin 2.5 --split-signal-samples
plot_mc_mc.py -a dplus_spg.yaml -v "Dmeson_m" --nology
```

### BKG comparison

Makes plots with 'Loose Inclusive' background templates and compares them to the SR-like background templates.

```
generate_wplusd_config.py -a dplus_bkg -s bkg_comparison --extra-rebin 2.5 --split-signal-samples
plot_mc_mc.py -a dplus_bkg.yaml -v "Dmeson_m" --nology
```

### data / MC plots without SPG and Loose Inclusive templates

Add the option `--replacement-samples` to replace some of the signal and background MC templates with SPG and Loose Inclusive templates.

```
generate_wplusd_config.py -a wplusd -s truth --extra-rebin 2.5 --split-signal-samples --replacement-samples
plot_data_mc.py -a wplusd.yaml -v Dmeson_m --nology
```

### Differential bins

Two additional options are used with the `generate_wplusd_config.py` script for differential bins:
- `--differential-bins`: Create additional plots for the differential pT(D) bins,
- `--truth-differential-bins`: Split the signal sample into five samples depending on the truth pT(D).
