#!/usr/bin/env python
from array import array
from charmplot.common import utils
from charmplot.control import channel
from charmplot.control import variable
import os
import ROOT


# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasLabels.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasUtils.C"))
ROOT.SetAtlasStyle()

# lumi
LUMI_RUN2 = 138965.16

# text precision
ROOT.gStyle.SetPaintTextFormat(".2f")

# files
f = ROOT.TFile("/global/cscratch1/sd/mmuskinj/charmpp/v8/stdm13_matrix_v2/wplusd_fit_wplusd_comparison/histograms.root", "READ")

# truth
f_truth = ROOT.TFile("/global/cscratch1/sd/mmuskinj/charmpp/v8/stdm13_truth_v1/truth/wplusd_truth_analysis/histograms.root", "READ")

# channels
channels = [
    "OS-SS_Dplus",
    "OS-SS_el_minus_Dplus",
    "OS-SS_el_plus_Dplus",
    "OS-SS_mu_minus_Dplus",
    "OS-SS_mu_plus_Dplus",
]

# samples
samples = [
    "MG_Wjets",
    "Powheg_Wjets",
    "Sherpa_Wjets",
]

# colors
colors = {
    "MG_Wjets": ROOT.kBlack,
    "Powheg_Wjets": ROOT.kBlue,
    "Sherpa_Wjets": ROOT.kRed,
}

# bin shift
bin_shift = {
    "MG_Wjets": 1.00,
    "Powheg_Wjets": 1.05,
    "Sherpa_Wjets": 0.95,
}

# variables
truth_pt = variable.Variable("D_pt", **{
    "label": "p^{truth}_{T}(D)",
    "unit": "GeV"})

# out folder
if not os.path.isdir("transfer_matrix"):
    os.mkdir("transfer_matrix")


def shifted_bins(xbins, r):
    out = array('d', [x * r for x in xbins])
    return out


# out file
f_out = ROOT.TFile("transfer_matrix/unfolding.root", "RECREATE")

# loop
for c in channels:

    # matrix
    h_tmp = {s: f.Get(f"{s}_emu_Matched_{c}_Dmeson_transfer_matrix") for s in samples}
    nbins = h_tmp["MG_Wjets"].GetNbinsX()
    xbins = array('d', [h_tmp["MG_Wjets"].GetXaxis().GetBinLowEdge(i) for i in range(1, nbins + 3)])
    h = {s: ROOT.TH2D(f"{s}_{c}_transfer_matrix", f"{s}_{c}_transfer_matrix", nbins + 1, xbins, nbins + 1, xbins) for s in samples}

    # differential bins
    h_pt_tmp = {s: f.Get(f"{s}_emu_Matched_{c}_Dmeson_differential_pt") for s in samples}
    h_pt = {s: ROOT.TH1D(f"{s}_{c}_differential_pt", f"{s}_{c}_differential_pt", nbins + 1, xbins) for s in samples}

    # truth differential
    h_pt_truth_tmp = {s: f_truth.Get(f"{s}_{c}_D_differential_pt") for s in samples}
    h_pt_truth = {s: ROOT.TH1D(f"{s}_{c}_truth_differential_pt", f"{s}_{c}_truth_differential_pt", nbins + 1, xbins) for s in samples}

    # calculate fiducial efficiency
    truth_projection_tmp = {s: h_tmp[s].ProjectionY(f"{s}_{c}_truth_projection_tmp", 0, nbins + 1) for s in samples}
    _ = [truth_projection_tmp[s].Scale(1. / LUMI_RUN2) for s in samples]
    truth_projection = {s: ROOT.TH1D(f"{s}_{c}_truth_projection", f"{s}_{c}_truth_projection", nbins + 1, xbins) for s in samples}

    # fill
    proxy_axis = {}
    fid_eff = {}
    fid_eff_gr = {}
    for s in samples:
        for i in range(1, nbins + 2):
            h_pt[s].SetBinContent(i, h_pt_tmp[s].GetBinContent(i))
            h_pt[s].SetBinError(i, h_pt_tmp[s].GetBinError(i))
            h_pt_truth[s].SetBinContent(i, h_pt_truth_tmp[s].GetBinContent(i))
            h_pt_truth[s].SetBinError(i, h_pt_truth_tmp[s].GetBinError(i))
            truth_projection[s].SetBinContent(i, truth_projection_tmp[s].GetBinContent(i))
            truth_projection[s].SetBinError(i, truth_projection_tmp[s].GetBinError(i))
            for j in range(1, nbins + 2):
                h[s].SetBinContent(i, j, h_tmp[s].GetBinContent(i, j))
                h[s].SetBinError(i, j, h_tmp[s].GetBinError(i, j))

        # axis title
        h[s].GetXaxis().SetTitle("p_{T}^{reco}(D) [GeV]")
        h[s].GetYaxis().SetTitle("p_{T}^{truth}(D) [GeV]")
        h[s].SetMarkerSize(1.4)

        h_pt[s].GetXaxis().SetTitle("p_{T}^{reco}(D) [GeV]")
        h_pt[s].GetYaxis().SetTitle("Entries / (bin)")
        h_pt[s].SetMarkerSize(1.4)
        h_pt[s].SetMarkerColor(colors[s])
        h_pt[s].SetLineColor(colors[s])

        h_pt_truth[s].GetXaxis().SetTitle("p_{T}^{truth}(D) [GeV]")
        h_pt_truth[s].GetYaxis().SetTitle("d#sigma / dp_{T}(D) [pb / bin]")
        h_pt_truth[s].SetMarkerSize(1.4)
        h_pt_truth[s].SetMarkerColor(colors[s])
        h_pt_truth[s].SetLineColor(colors[s])

        proxy_axis[s] = ROOT.TH1D(f"{s}_{c}_proxy_axis", f"{s}_{c}_proxy_axis", nbins + 1, shifted_bins(xbins, bin_shift[s]))
        proxy_axis[s].GetXaxis().SetMoreLogLabels()
        proxy_axis[s].SetLineWidth(0)
        proxy_axis[s].SetLineWidth(0)
        proxy_axis[s].SetMinimum(0.00)
        proxy_axis[s].SetMaximum(0.025)
        proxy_axis[s].SetMarkerSize(1.4)
        proxy_axis[s].SetMarkerColor(colors[s])

        # calculate fiducial efficiency
        fid_eff[s] = ROOT.TEfficiency(truth_projection[s], h_pt_truth[s])
        fid_eff_gr[s] = fid_eff[s].CreateGraph()
        fid_eff_gr[s].SetMarkerColor(colors[s])
        fid_eff_gr[s].SetLineColor(colors[s])

    f_out.cd()
    for s in samples:
        h[s].Write()
        h_pt[s].Write()
        h_pt_truth[s].Write()
        fid_eff_gr[s].Write(f"{s}_{c}_fid_eff")

    # -------------------
    # draw matrix
    # -------------------
    canv1 = ROOT.TCanvas(f"{c}_matrix", f"{c}_matrix", 1000, 800)
    canv1.SetRightMargin(0.15)
    canv1.SetLogy()
    canv1.SetLogx()
    h["MG_Wjets"].GetXaxis().SetMoreLogLabels()
    h["MG_Wjets"].GetYaxis().SetMoreLogLabels()
    h["MG_Wjets"].Draw("text colz error")

    # ATLAS label
    ROOT.ATLASLabel(0.18, 0.90, "Internal", 1)
    ROOT.myText(0.18, 0.84, 1, "#sqrt{s} = 13 TeV")
    ROOT.myText(0.18, 0.78, 1, "139 fb^{-1}")
    ROOT.myText(0.46, 0.24, 1, "W#rightarrowl#nu+D, D#rightarrowK#pi#pi")
    ROOT.myText(0.46, 0.18, 1, c.replace("OS-SS_", ""))

    # save
    ROOT.gPad.RedrawAxis()
    canv1.Print(f"transfer_matrix/{c}_matrix.pdf")

    # normalize
    _ = [h[s].Scale(100. / h[s].GetSumOfWeights()) for s in samples]
    canv1.Print(f"transfer_matrix/{c}_matrix_NORM.pdf")

    # legend
    leg = ROOT.TLegend(0.65, 0.90 - 0.05 * len(samples), 0.90, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(28)
    leg.SetTextFont(43)
    _ = [leg.AddEntry(fid_eff_gr[s], s, "pe") for s in samples]

    # channel ojbect
    chan = channel.Channel(c, ["W#rightarrowl#nu+D, D#rightarrowK#pi#pi", c], "", [], [])

    # -------------------
    # draw reco
    # -------------------
    canv2 = ROOT.TCanvas(f"{c}_reco", f"{c}_reco", 1000, 800)
    canv2.SetLogx()
    h_pt["MG_Wjets"].SetMaximum(2 * h_pt["MG_Wjets"].GetMaximum())
    h_pt["MG_Wjets"].GetXaxis().SetMoreLogLabels()
    h_pt["MG_Wjets"].Draw("text hist")
    for s in samples[1:]:
        h_pt[s].Draw("hist same")
    leg.Draw()

    # ATLAS label
    ROOT.ATLASLabel(0.18, 0.90, "Internal", 1)
    ROOT.myText(0.18, 0.84, 1, "#sqrt{s} = 13 TeV, 139 fb^{-1}")
    ROOT.myText(0.18, 0.78, 1, "W#rightarrowl#nu+D, D#rightarrowK#pi#pi")
    ROOT.myText(0.18, 0.72, 1, c)

    # save
    ROOT.gPad.RedrawAxis()
    canv2.Print(f"transfer_matrix/{c}_reco.pdf")

    # -------------------
    # draw truth
    # -------------------
    canv3 = ROOT.TCanvas(f"{c}_truth", f"{c}_truth", 1000, 800)
    canv3.SetLogx()
    h_pt_truth["MG_Wjets"].SetMaximum(2 * h_pt_truth["MG_Wjets"].GetMaximum())
    h_pt_truth["MG_Wjets"].GetXaxis().SetMoreLogLabels()
    h_pt_truth["MG_Wjets"].Draw("text hist")
    for s in samples[1:]:
        h_pt_truth[s].Draw("hist same")
    leg.Draw()

    # ATLAS label
    ROOT.ATLASLabel(0.18, 0.90, "Internal", 1)
    ROOT.myText(0.18, 0.84, 1, "#sqrt{s} = 13 TeV")
    ROOT.myText(0.18, 0.78, 1, "W#rightarrowl#nu+D, D#rightarrowK#pi#pi")
    ROOT.myText(0.18, 0.72, 1, c)

    # save
    ROOT.gPad.RedrawAxis()
    canv3.Print(f"transfer_matrix/{c}_truth.pdf")

    # -------------------
    # draw fiducial efficiency
    # -------------------
    ROOT.gStyle.SetPaintTextFormat(".5f")
    canv4 = utils.make_canvas_mc_ratio(proxy_axis[samples[0]], truth_pt, chan, "Ratio", x=800, y=800, events="fiducial efficiency")
    canv4.proxy_dn.GetXaxis().SetLabelOffset(-0.03)
    canv4.pad1.cd()
    canv4.pad1.SetLogx()
    for s in samples:
        for i in range(1, nbins + 2):
            proxy_axis[s].SetBinContent(i, fid_eff_gr[s].GetY()[i - 1])
            proxy_axis[s].SetBinError(i, 0)
            proxy_axis[s].Draw("hist text same")
    _ = [fid_eff_gr[s].Draw("pe") for s in reversed(samples)]
    leg.Draw()

    # ratio
    canv4.pad2.cd()
    canv4.pad2.SetLogx()
    canv4.set_ratio_range(0.91, 1.09, override=True)
    fid_eff_gr_ratio = {s: fid_eff_gr[s].Clone(f"{fid_eff_gr[s].GetName()}_ratio") for s in samples}
    for s in samples:
        for i in range(fid_eff_gr_ratio[s].GetN()):
            fid_eff_gr_ratio[s].GetY()[i] /= fid_eff_gr[samples[0]].GetY()[i]
            fid_eff_gr_ratio[s].GetEYhigh()[i] /= fid_eff_gr[samples[0]].GetY()[i]
            fid_eff_gr_ratio[s].GetEYlow()[i] /= fid_eff_gr[samples[0]].GetY()[i]
        _ = [fid_eff_gr_ratio[s].Draw("pe") for s in reversed(samples)]

    # save
    ROOT.gPad.RedrawAxis()
    canv4.print(f"transfer_matrix/{c}_fid_eff.pdf")

# close file
f_out.Close()
