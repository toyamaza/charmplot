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

FOLDER = "/global/cscratch1/sd/mmuskinj/charmpp/charm_frag_v1/"

LUMI = 58450.1


def main():

    f_reco = ROOT.TFile(os.path.join(FOLDER, "reco_v1/wplusd_charm_frag/histograms.root"))
    f_truth = ROOT.TFile(os.path.join(FOLDER, "truth_v1/wplusd_truth_charm_frag/histograms.root"))

    h_reco = f_reco.Get("Sherpa2211_Wjets_OS-SS_0tag_Dplus_Dmeson_transfer_matrix_zt")
    h_truth = f_truth.Get("Sherpa2211_Wjets_OS-SS_0tag_Dplus_D_jet_zt")

    for i in range(1, h_truth.GetNbinsX() + 1):
        for j in range(1, h_reco.GetNbinsX() + 1):
            h_reco.SetBinContent(i, j, h_reco.GetBinContent(i, j) / (LUMI * h_truth.GetBinContent(i)))

    # normalized matrix
    integral = h_reco.GetSumOfWeights()
    for i in range(1, h_reco.GetNbinsX() + 1):
        for j in range(1, h_reco.GetNbinsX() + 1):
            z = h_reco.GetBinContent(i, j)
            h_reco.SetBinContent(i, j, 100 * z / integral)

    ROOT.gStyle.SetPaintTextFormat("1.1e")
    ROOT.gStyle.SetPalette(ROOT.kCherry)
    ROOT.TColor.InvertPalette()
    canv = ROOT.TCanvas("response_matrix", "response_matrix", 1000, 800)
    canv.SetRightMargin(0.18)
    canv.SetLeftMargin(0.18)
    canv.SetBottomMargin(0.20)
    h_reco.GetZaxis().SetTitle("Normalized Response Matrix [%]")
    h_reco.GetXaxis().SetTitle("z_{T}^{D/jet} (reco)")
    h_reco.GetYaxis().SetTitle("z_{T}^{D/jet} (truth)")
    h_reco.GetZaxis().SetTitleSize(h_reco.GetZaxis().GetTitleSize() * 0.9)
    h_reco.GetZaxis().SetLabelSize(h_reco.GetZaxis().GetLabelSize() * 0.97)
    h_reco.SetMinimum(1e-2)
    h_reco.Draw("colz")

    # ATLAS label
    ROOT.ATLASLabel(0.20, 0.90, "Simulation Internal", 1, 0.04)
    ROOT.myText(0.20, 0.85, 1, "#sqrt{s} = 13 TeV", 0.04)
    ROOT.myText(0.20, 0.80, 1, "W+D(#rightarrowK#pi#pi)", 0.04)

    # save
    ROOT.gPad.RedrawAxis()
    canv.Print("transfer_matrix.pdf")


if __name__ == "__main__":

    # run
    main()
