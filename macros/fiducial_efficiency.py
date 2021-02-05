#!/usr/bin/env python
from ctypes import c_double
import os
import ROOT

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.SetAtlasStyle()

# lumi
lumi = 138965.16

# files
path = "/global/cscratch1/sd/mmuskinj/charmpp/v8/fid_eff"
f_LO_MG = ROOT.TFile(f"{path}/MG_Wjets_emu.root", "OPEN")
f_MG_FxFx = ROOT.TFile(f"{path}/MG_FxFx_Wjets_emu.root", "OPEN")
f_Sherpa2210 = ROOT.TFile(f"{path}/Sh_2210_Wjets.root", "OPEN")
f_Powheg = ROOT.TFile(f"{path}/Ph_Wjets_emu.root", "OPEN")

# muon channel
mu_name = "mu_SR_Dplus_OS_MatchedFid__Dmeson_m"

# OS - SS
h_LO_MG = f_LO_MG.Get(mu_name).Clone(f"{f_LO_MG}_LO_MG")
h_MG_FxFx = f_MG_FxFx.Get(mu_name).Clone(f"{f_LO_MG}_LO_MG")
h_Sherpa2210 = f_Sherpa2210.Get(mu_name).Clone(f"{f_LO_MG}_LO_MG")
h_Powheg = f_Powheg.Get(mu_name).Clone(f"{f_LO_MG}_LO_MG")
h_LO_MG.Add(f_LO_MG.Get(mu_name.replace("OS", "SS")), -1)
h_MG_FxFx.Add(f_MG_FxFx.Get(mu_name.replace("OS", "SS")), -1)
h_Sherpa2210.Add(f_Sherpa2210.Get(mu_name.replace("OS", "SS")), -1)
h_Powheg.Add(f_Powheg.Get(mu_name.replace("OS", "SS")), -1)

# integrals
err_LO_MG = c_double()
err_MG_FxFx = c_double()
err_Sherpa2210 = c_double()
err_Powheg = c_double()
int_LO_MG = h_LO_MG.IntegralAndError(0, h_LO_MG.GetNbinsX() + 2, err_LO_MG)
int_MG_FxFx = h_MG_FxFx.IntegralAndError(0, h_MG_FxFx.GetNbinsX() + 2, err_MG_FxFx)
int_Sherpa2210 = h_Sherpa2210.IntegralAndError(0, h_Sherpa2210.GetNbinsX() + 2, err_Sherpa2210)
int_Powheg = h_Powheg.IntegralAndError(0, h_Powheg.GetNbinsX() + 2, err_Powheg)

# print
print(f"LO_MG: {int_LO_MG/lumi} {err_LO_MG.value/lumi}")
print(f"Sherpa2210: {int_Sherpa2210/lumi} {err_Sherpa2210.value/lumi}")
print(f"MG_FxFx: {int_MG_FxFx/lumi} {err_MG_FxFx.value/lumi}")
print(f"Powheg: {int_Powheg/lumi} {err_Powheg.value/lumi}")
