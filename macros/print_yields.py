#!/usr/bin/env python
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
f = ROOT.TFile("/global/cscratch1/sd/mmuskinj/charmpp/v10_Sherpa2211_fit/fit_6/wplusd_nosys/histograms.root", "READ")

# for loop
for i in range(1, 6):
    string = ""
    for lep in ["lep_minus", "lep_plus"]:
        for charge in ["OS", "SS"]:
            for btag in ["0tag", "1tag"]:
                hname = f"Data_{charge}_{lep}_{btag}_Dplus_pt_bin{i}_Dmeson_m"
                h = f.Get(hname)
                string += f"\t{h.Integral()}"
    print(string)
