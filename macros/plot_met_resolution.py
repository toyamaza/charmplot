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

# input
f = ROOT.TFile("wplusd_comparison/histograms.root", "READ")
h = f.Get("Sherpa2211_WplusD_OS-SS_0tag_Dplus_Matched_truth_met_met_resolution")

# make plot
canv = ROOT.TCanvas("met_reso", "met_res", 800, 600)
canv.SetRightMargin(0.20)
h.Draw("COLZ")
h.GetXaxis().SetTitle("MET^{reco} [GeV]")
h.GetYaxis().SetTitle("MET^{truth} [GeV]")
h.GetZaxis().SetTitle("Entries")
h.GetZaxis().SetTitleOffset(1.2 * h.GetZaxis().GetTitleOffset())

# ATLAS label
ROOT.ATLASLabel(0.18, 0.90, "Internal", 0)
ROOT.myText(0.18, 0.84, 0, "#sqrt{s} = 13 TeV, 139 fb^{-1}")
ROOT.myText(0.18, 0.78, 0, "OS-SS, W+D signal")

# lines
line1 = ROOT.TLine(30, 0, 30, 70)
line1.SetLineColor(0)
line1.SetLineStyle(2)
line2 = ROOT.TLine(0, 30, 100, 30)
line2.SetLineColor(0)
line2.SetLineStyle(2)
line1.Draw()
line2.Draw()

canv.Print("met_resolution.pdf")
