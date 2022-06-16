#!/usr/bin/env python
from ctypes import c_double
import os
import ROOT


# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasLabels.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasUtils.C"))
ROOT.SetAtlasStyle()

# input file
f = ROOT.TFile("truth/wplusd_truth_analysis/histograms.root", "READ")

# channel
for s in ["Sherpa2211_Wjets", "Sherpa2211_WplusD", "MGFxFx_Wjets", "MGFxFx_WplusD", "MGPy8EG_NLO_WplusD", "MG_Wjets"]:
    for d in ["Dplus", "Dstar"]:
        print(f"{s} {d}")
        for c in ["OS-SS_el_minus", "OS-SS_el_plus", "OS-SS_mu_minus", "OS-SS_mu_plus", "OS-SS"]:
            h = f.Get(f"{s}_{c}_{d}_D_pt_truth")
            err = c_double(0)
            integral = h.IntegralAndError(0, h.GetNbinsX() + 1, err)
            if d == "Dstar":
                err.value = err.value / 0.677
                integral = integral / 0.677
            print(f"{integral}\t{err.value}")
