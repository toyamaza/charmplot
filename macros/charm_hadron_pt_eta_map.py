#!/usr/bin/env python
import os
import ROOT

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.SetAtlasStyle()

# input file
f = ROOT.TFile("/global/cscratch1/sd/mmuskinj/charmpp/v8/truth/truth/charm_truth/histograms.root", "READ")

h_Dplus_plus = f.Get("WPlusD_OS_minus_SS_el_minus_Dplus_D_pt_eta").Clone("Dplus_plus")
h_Dplus_plus.Add(f.Get("WPlusD_OS_minus_SS_mu_minus_Dplus_D_pt_eta"))

h_Dplus_minus = f.Get("WPlusD_OS_minus_SS_el_plus_Dplus_D_pt_eta").Clone("Dplus_minus")
h_Dplus_minus.Add(f.Get("WPlusD_OS_minus_SS_mu_plus_Dplus_D_pt_eta"))

h_Dstar_plus = f.Get("WPlusD_OS_minus_SS_el_minus_Dstar_D_pt_eta").Clone("Dstar_plus")
h_Dstar_plus.Add(f.Get("WPlusD_OS_minus_SS_mu_minus_Dstar_D_pt_eta"))

h_Dstar_minus = f.Get("WPlusD_OS_minus_SS_el_plus_Dstar_D_pt_eta").Clone("Dstar_minus")
h_Dstar_minus.Add(f.Get("WPlusD_OS_minus_SS_mu_plus_Dstar_D_pt_eta"))

f_out = ROOT.TFile("CharmMaps.root", "RECREATE")
f_out.cd()
h_Dplus_plus.Write()
h_Dplus_minus.Write()
h_Dstar_plus.Write()
h_Dstar_minus.Write()

h_Dplus_plus.ProjectionX(f"{h_Dplus_plus.GetName()}_projX").Write()
h_Dplus_minus.ProjectionX(f"{h_Dplus_minus.GetName()}_projX").Write()
h_Dstar_plus.ProjectionX(f"{h_Dstar_plus.GetName()}_projX").Write()
h_Dstar_minus.ProjectionX(f"{h_Dstar_minus.GetName()}_projX").Write()
