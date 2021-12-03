#!/usr/bin/env python
from array import array

from charmplot.common import utils
from charmplot.control import channel
from charmplot.control import variable
import numpy as np
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

# colors
colors = {
    "MG_Wjets": ROOT.kBlack,
    "Powheg_Wjets": ROOT.kBlue,
    "Sherpa_Wjets": ROOT.kRed,
    "Sherpa2211_Wjets": ROOT.kGreen,
    "MGFxFx_Wjets": ROOT.kViolet,
    "Sherpa2211_WplusD": ROOT.kBlue,
    "MGPy8EG_NLO_WplusD": ROOT.kRed,
}

# bin shift
bin_shift = {
    "MG_Wjets": 1.00,
    "Powheg_Wjets": 1.05,
    "Sherpa_Wjets": 0.95,
    "Sherpa2211_Wjets": 0.925,
    "MGFxFx_Wjets": 1.25,
    "Sherpa2211_WplusD": 1.00,
    "MGPy8EG_NLO_WplusD": 1.00,
}

# legend
legend_names = {
    "MG_Wjets": "LO MadGraph",
    "Powheg_Wjets": "Powheg",
    "Sherpa_Wjets": "Sherpa2.2.1",
    "Sherpa2211_Wjets": "Sh.2.11 W+jets",
    "MGFxFx_Wjets": "MG FxFx",
    "Sherpa2211_WplusD": "Sh.2.11 W+D",
    "MGPy8EG_NLO_WplusD": "NLO MG W+D",
}

# variables
truth_pt = variable.Variable("D_pt", **{
    "label": "p^{truth}_{T}(D)",
    "unit": "GeV"})

reco_pt = variable.Variable("D_pt", **{
    "label": "p^{reco}_{T}(D)",
    "unit": "GeV"})


def shifted_bins(xbins, r):
    out = array('d', [x * r for x in xbins])
    return out


def main(options, args):

    # out folder
    if not os.path.isdir(options.output):
        os.mkdir(options.output)

    # out file
    f_out = ROOT.TFile(f"{options.output}/unfolding.root", "RECREATE")

    # input files
    f = ROOT.TFile(options.reco, "READ")
    f_truth = ROOT.TFile(options.truth, "READ")

    # samples
    samples = options.samples.split(",")
    print(f"Got samples: {samples}")

    # channels
    channels = options.channels.split(",")
    print(f"Got channels: {samples}")

    # systematics
    systematics = [""]
    if options.sherpa_pdf:
        systematics += [f"GEN_ME_ONLY_MUR1_MUF1_PDF{N}" for N in range(303200, 303301)]
    elif options.sherpa_qcd:
        systematics += [
            "GEN_ME_ONLY_MUR05_MUF05_PDF303200_PSMUR05_PSMUF05",
            "GEN_ME_ONLY_MUR05_MUF1_PDF303200_PSMUR05_PSMUF1",
            "GEN_ME_ONLY_MUR1_MUF05_PDF303200_PSMUR1_PSMUF05",
            "GEN_ME_ONLY_MUR1_MUF2_PDF303200_PSMUR1_PSMUF2",
            "GEN_ME_ONLY_MUR2_MUF1_PDF303200_PSMUR2_PSMUF1",
            "GEN_ME_ONLY_MUR2_MUF2_PDF303200_PSMUR2_PSMUF2",
        ]
    elif options.sherpa_as:
        systematics += [
            "GEN_ME_ONLY_MUR1_MUF1_PDF270000",
            "GEN_ME_ONLY_MUR1_MUF1_PDF269000",
        ]
    elif options.sherpa_ew:
        systematics += [
            "GEN_ME_ONLY_MUR1_MUF1_PDF303200_ASSEW",
            "GEN_ME_ONLY_MUR1_MUF1_PDF303200_ASSEWLO1",
            "GEN_ME_ONLY_MUR1_MUF1_PDF303200_ASSEWLO1LO2",
            "GEN_ME_ONLY_MUR1_MUF1_PDF303200_ASSEWLO1LO2LO3",
            "GEN_ME_ONLY_MUR1_MUF1_PDF303200_EXPASSEW",
            "GEN_ME_ONLY_MUR1_MUF1_PDF303200_EXPASSEWLO1",
            "GEN_ME_ONLY_MUR1_MUF1_PDF303200_EXPASSEWLO1LO2",
            "GEN_ME_ONLY_MUR1_MUF1_PDF303200_EXPASSEWLO1LO2LO3",
            "GEN_ME_ONLY_MUR1_MUF1_PDF303200_MULTIASSEW",
            "GEN_ME_ONLY_MUR1_MUF1_PDF303200_MULTIASSEWLO1",
            "GEN_ME_ONLY_MUR1_MUF1_PDF303200_MULTIASSEWLO1LO2",
            "GEN_ME_ONLY_MUR1_MUF1_PDF303200_MULTIASSEWLO1LO2LO3",
        ]

    # samples sys
    samples_sys = []
    for s in samples:
        for syst in systematics:
            if syst:
                if not s == "Sherpa2211_WplusD":
                    continue
                syst = f"_{syst}"
            samples_sys += [f"{s}{syst}"]
    print(f"Samples sys: {samples_sys}")

    # loop
    for c in channels:

        print(f"Doing channel {c}")

        # reco
        h_tmp = {}
        h = {}
        h_pt_tmp = {}
        h_pt = {}

        # truth
        h_pt_truth_tmp = {}
        h_pt_truth = {}
        truth_projection_tmp = {}
        truth_projection = {}
        h_fid_eff = {}
        h_fid_eff_inv = {}

        for s in samples_sys:

            # systematic and special cases
            sys = s.split("_GEN")

            # special name for truth-matched reco samples
            if "WplusD" not in s and "Sherpa2211" not in s:
                name_reco = sys[0] + "_emu_Matched"
            else:
                name_reco = sys[0]
            name = sys[0]

            # systematics
            if len(sys) > 1:
                sys = "_GEN" + sys[-1]
            else:
                sys = ""

            # histogram name
            print(f"{name_reco}_{c.replace('_Kpipi', '')}_Dmeson_transfer_matrix{sys}")

            h_tmp[s] = f.Get(f"{name_reco}_{c.replace('_Kpipi', '')}_Dmeson_transfer_matrix{sys}")
            nbins = h_tmp[s].GetNbinsX()
            xbins = array('d', [h_tmp[s].GetXaxis().GetBinLowEdge(i) for i in range(1, nbins + 3)])
            xbins[-1] = xbins[-2] + 2 * xbins[-3]
            h[s] = ROOT.TH2D(f"{s}_{c}_transfer_matrix", f"{s}_{c}_transfer_matrix", nbins + 1, xbins, nbins + 1, xbins)

            # differential bins
            h_pt_tmp[s] = f.Get(f"{name_reco}_{c.replace('_Kpipi', '')}_Dmeson_differential_pt{sys}")
            h_pt[s] = ROOT.TH1D(f"{s}_{c}_differential_pt", f"{s}_{c}_differential_pt", nbins + 1, xbins)

            # truth histogram name
            print(f"{name}_{c.replace('_0tag', '')}_D_differential_pt{sys}")

            h_pt_truth_tmp[s] = f_truth.Get(f"{name}_{c.replace('_0tag', '')}_D_differential_pt{sys}")
            h_pt_truth[s] = ROOT.TH1D(f"{s}_{c}_truth_differential_pt", f"{s}_{c}_truth_differential_pt", nbins + 1, xbins)

            # calculate fiducial efficiency
            truth_projection_tmp[s] = h_tmp[s].ProjectionY(f"{s}_{c}_truth_projection_tmp", 0, nbins + 1)
            truth_projection_tmp[s].Scale(1. / LUMI_RUN2)
            truth_projection[s] = ROOT.TH1D(f"{s}_{c}_truth_projection", f"{s}_{c}_truth_projection", nbins + 1, xbins)

            # # account for forced decay samples
            # if "Kpipi" in c:
            #     for s in h_pt_truth:
            #         if "WplusD" in s:
            #             h_pt_truth_tmp[s].Scale(0.096)

            # calculate fiducial efficiency per bin
            h_fid_eff[s] = ROOT.TH2D(f"{s}_{c}_fid_eff", f"{s}_{c}_fid_eff", nbins + 1, xbins, nbins + 1, xbins)
            h_fid_eff_inv[s] = ROOT.TH2D(f"{s}_{c}_fid_eff_inv", f"{s}_{c}_fid_eff_inv", nbins + 1, xbins, nbins + 1, xbins)

        # numpy matrix
        np_matrix = np.identity(nbins + 1)
        np_truth = np.ones(nbins + 1)

        # fill
        proxy_axis = {}
        fid_eff = {}
        fid_eff_gr = {}
        fid_eff_inclusive = {}
        for s in samples_sys:
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

            h_pt_truth[s].GetXaxis().SetTitle("p_{T}^{truth}(D) [GeV]")
            h_pt_truth[s].GetYaxis().SetTitle("d#sigma / dp_{T}(D) [pb / bin]")
            h_pt_truth[s].SetMarkerSize(1.4)

            h_fid_eff[s].GetXaxis().SetTitle("p_{T}^{reco}(D) [GeV]")
            h_fid_eff[s].GetYaxis().SetTitle("p_{T}^{truth}(D) [GeV]")
            h_fid_eff[s].SetMarkerSize(1.4)
            h_fid_eff[s].SetMarkerColor(ROOT.kWhite)

            h_fid_eff_inv[s].GetXaxis().SetTitle("p_{T}^{reco}(D) [GeV]")
            h_fid_eff_inv[s].GetYaxis().SetTitle("p_{T}^{truth}(D) [GeV]")
            h_fid_eff_inv[s].SetMarkerSize(1.4)
            h_fid_eff_inv[s].SetMarkerColor(ROOT.kWhite)

            # calculate fiducial efficiency
            fid_eff[s] = ROOT.TEfficiency(truth_projection[s], h_pt_truth[s])
            fid_eff_gr[s] = fid_eff[s].CreateGraph()

            # color
            if "_GEN_" not in s:
                h_pt[s].SetMarkerColor(colors[s])
                h_pt[s].SetLineColor(colors[s])
                h_pt_truth[s].SetMarkerColor(colors[s])
                h_pt_truth[s].SetLineColor(colors[s])
                fid_eff_gr[s].SetMarkerColor(colors[s])
                fid_eff_gr[s].SetLineColor(colors[s])

                # proxy axis for fiducil efficinecy
                proxy_axis[s] = ROOT.TH1D(f"{s}_{c}_proxy_axis", f"{s}_{c}_proxy_axis", nbins + 1, shifted_bins(xbins, bin_shift[s]))
                proxy_axis[s].GetXaxis().SetNoExponent()
                proxy_axis[s].SetLineWidth(0)
                proxy_axis[s].SetLineWidth(0)
                proxy_axis[s].SetMinimum(0.00)
                proxy_axis[s].SetMaximum(0.025)
                if "Kpipi" in c:
                    proxy_axis[s].SetMaximum(0.025 / 0.092)
                proxy_axis[s].SetMarkerSize(1.4)
                proxy_axis[s].SetMarkerColor(colors[s])
            else:
                h_pt[s].SetMarkerColor(colors[s.split("_GEN_")[0]] + 2)
                h_pt[s].SetLineColor(colors[s.split("_GEN_")[0]] + 2)
                h_pt_truth[s].SetMarkerColor(colors[s.split("_GEN_")[0]] + 2)
                h_pt_truth[s].SetLineColor(colors[s.split("_GEN_")[0]] + 2)
                fid_eff_gr[s].SetMarkerColor(colors[s.split("_GEN_")[0]] + 2)
                fid_eff_gr[s].SetLineColor(colors[s.split("_GEN_")[0]] + 2)

            if "_GEN_" not in s:
                # calculate fiducial efficiency per bin
                for i in range(1, nbins + 2):
                    for j in range(1, nbins + 2):
                        tmp_num = ROOT.TH1D(f"{h[s].GetName()}_tmp_{i}_{j}", f"{h[s].GetName()}_tmp_{i}_{j}", 1, 0, 1)
                        tmp_num.SetBinContent(1, h[s].GetBinContent(i, j))
                        tmp_num.SetBinError(1, h[s].GetBinError(i, j))
                        tmp_den = ROOT.TH1D(f"{h[s].GetName()}_tmp_{j}", f"{h[s].GetName()}_tmp_{j}", 1, 0, 1)
                        tmp_den.SetBinContent(1, LUMI_RUN2 * h_pt_truth[s].GetBinContent(j))
                        tmp_den.SetBinError(1, LUMI_RUN2 * h_pt_truth[s].GetBinError(j))
                        tmp_eff = ROOT.TEfficiency(tmp_num, tmp_den)
                        h_fid_eff[s].SetBinContent(i, j, 100 * tmp_eff.GetEfficiency(1))
                        np_matrix[i - 1][j - 1] = tmp_eff.GetEfficiency(1)
                        np_truth[j - 1] = LUMI_RUN2 * h_pt_truth[s].GetBinContent(j)
                        if tmp_eff.GetEfficiency(1) > 0:
                            h_fid_eff[s].SetBinError(i, j, 100 * ((tmp_eff.GetEfficiencyErrorUp(1) +
                                                                   tmp_eff.GetEfficiencyErrorLow(1)) / 2) / tmp_eff.GetEfficiency(1))
                        print(f"{i} {j} {tmp_eff.GetEfficiency(1)} {(tmp_eff.GetEfficiencyErrorUp(1) + tmp_eff.GetEfficiencyErrorLow(1)) / 2}")

            # inclusive efficiency
            inclusive_num = truth_projection[s].Clone(f"{truth_projection[s].GetName()}_inclusive")
            inclusive_den = h_pt_truth[s].Clone(f"{h_pt_truth[s].GetName()}_inclusive")
            inclusive_num.Rebin(inclusive_num.GetNbinsX())
            inclusive_den.Rebin(inclusive_den.GetNbinsX())
            fid_eff_inclusive[s] = ROOT.TEfficiency(inclusive_num, inclusive_den)

        f_out.cd()

        # -------------------
        # draw matrix
        # -------------------
        canv1 = ROOT.TCanvas(f"{c}_matrix", f"{c}_matrix", 1000, 800)
        canv1.SetRightMargin(0.15)
        canv1.SetLogy()
        canv1.SetLogx()
        h[samples[0]].GetXaxis().SetMoreLogLabels()
        h[samples[0]].GetXaxis().SetNoExponent()
        h[samples[0]].GetYaxis().SetMoreLogLabels()
        h[samples[0]].Draw("text colz error")

        # ATLAS label
        ROOT.ATLASLabel(0.18, 0.90, "Internal", 1)
        ROOT.myText(0.18, 0.84, 1, "#sqrt{s} = 13 TeV")
        ROOT.myText(0.18, 0.78, 1, "139 fb^{-1}")
        if "Dplus" in c:
            ROOT.myText(0.50, 0.24, 1, "W#rightarrowl#nu+D, D#rightarrowK#pi#pi")
        elif "Dstar" in c:
            ROOT.myText(0.50, 0.24, 1, "W#rightarrowl#nu+D*, D*#rightarrow(K#pi)#pi")
        ROOT.myText(0.50, 0.18, 1, c.replace("OS-SS_", ""))

        # save
        ROOT.gPad.RedrawAxis()
        canv1.Print(f"{options.output}/{c}_matrix.pdf")

        # normalize
        _ = [h[s].Scale(100. / h[s].GetSumOfWeights()) for s in samples_sys]
        canv1.Print(f"{options.output}/{c}_matrix_NORM.pdf")

        # legend
        leg = ROOT.TLegend(0.65, 0.90 - 0.05 * len(samples), 0.90, 0.90)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(28)
        leg.SetTextFont(43)
        _ = [leg.AddEntry(fid_eff_gr[s], legend_names[s], "pe") for s in samples_sys if "_GEN_" not in s]

        # channel ojbect
        if "Dplus" in c:
            chan = channel.Channel(c, ["W#rightarrowl#nu+D, D#rightarrowK#pi#pi", c], "", [], [])
        elif "Dstar" in c:
            chan = channel.Channel(c, ["W#rightarrowl#nu+D*, D*#rightarrow(K#pi)#pi", c], "", [], [])

        # -------------------
        # draw reco
        # -------------------
        if options.sherpa_pdf:
            sherpa_nominal = None
            sherpa_sys = []
            for key, h_sys in h_pt.items():
                if "Sherpa2211_WplusD" not in key:
                    continue
                if key == "Sherpa2211_WplusD":
                    sherpa_nominal = h_sys
                else:
                    sherpa_sys += [h_sys]
            reco_sys_band, reco_sys_band_ratio = utils.make_pdf_err(sherpa_nominal, sherpa_sys, "NNPDF30_nnlo_as_0118_hessian")
        h_pt[samples_sys[0]].GetXaxis().SetNoExponent()

        canv2 = utils.make_canvas_mc_ratio(h_pt[samples_sys[0]], reco_pt, chan, "Ratio", x=800, y=800, events="d#sigma / dp_{T}(D) [pb / bin]", suffix="reco")
        canv2.pad1.cd()
        canv2.pad1.SetLogx()
        if options.sherpa_pdf:
            reco_sys_band.SetFillColor(colors["Sherpa2211_WplusD"] + 2)
            reco_sys_band.Draw("e2")
        for s in reversed(samples_sys):
            if (len(samples) > 1 or options.sherpa_pdf) and "_GEN_" in s:
                continue
            h_pt[s].Draw("hist same")
        leg.Draw()

        # set maximum after creating legend
        canv2.proxy_up.SetMaximum(2 * max([h_pt[s].GetMaximum() for s in h_pt]))

        # ratio
        ROOT.gPad.RedrawAxis()
        canv2.pad2.cd()
        canv2.pad2.SetLogx()
        if options.sherpa_pdf:
            reco_sys_band_ratio.SetFillColor(colors["Sherpa2211_WplusD"] + 2)
            reco_sys_band_ratio.Draw("e2")
        denum = h_pt[samples_sys[0]].Clone("reco_denum")
        for i in range(1, denum.GetNbinsX() + 1):
            denum.SetBinError(i, 0)
        ratios = []
        for s in reversed(samples_sys):
            if (len(samples) > 1 or options.sherpa_pdf) and "_GEN_" in s:
                continue
            h_pt_r = h_pt[s].Clone(f"{h_pt[s].GetName()}_r")
            h_pt_r.Divide(denum)
            h_pt_r.Draw("hist same")
            ratios += [h_pt_r]

        # save
        ROOT.gPad.RedrawAxis()
        canv2.set_ratio_range(0.71, 1.29, override=True)
        canv2.print(f"{options.output}/{c}_reco.pdf")
        h_pt[samples[0]].Write()

        # -------------------
        # draw truth
        # -------------------
        if options.sherpa_pdf:
            sherpa_nominal = None
            sherpa_sys = []
            for key, h_sys in h_pt_truth.items():
                if "Sherpa2211_WplusD" not in key:
                    continue
                if key == "Sherpa2211_WplusD":
                    sherpa_nominal = h_sys
                else:
                    sherpa_sys += [h_sys]
            truth_sys_band, truth_sys_band_ratio = utils.make_pdf_err(sherpa_nominal, sherpa_sys, "NNPDF30_nnlo_as_0118_hessian")
        h_pt_truth[samples_sys[0]].GetXaxis().SetNoExponent()

        canv3 = utils.make_canvas_mc_ratio(h_pt_truth[samples_sys[0]], truth_pt, chan, "Ratio", x=800,
                                           y=800, events="d#sigma / dp_{T}(D) [pb / bin]", suffix="truth")
        canv3.pad1.cd()
        canv3.pad1.SetLogx()
        if options.sherpa_pdf:
            truth_sys_band.SetFillColor(colors["Sherpa2211_WplusD"] + 2)
            truth_sys_band.Draw("e2")
        for s in reversed(samples_sys):
            if (len(samples) > 1 or options.sherpa_pdf) and "_GEN_" in s:
                continue
            h_pt_truth[s].Draw("hist same")
        leg.Draw()

        # set maximum after creating legend
        canv3.proxy_up.SetMaximum(2 * max([h_pt_truth[s].GetMaximum() for s in h_pt_truth]))

        # ratio
        ROOT.gPad.RedrawAxis()
        canv3.pad2.cd()
        canv3.pad2.SetLogx()
        if options.sherpa_pdf:
            truth_sys_band_ratio.SetFillColor(colors["Sherpa2211_WplusD"] + 2)
            truth_sys_band_ratio.Draw("e2")
        denum = h_pt_truth[samples_sys[0]].Clone("truth_denum")
        for i in range(1, denum.GetNbinsX() + 1):
            denum.SetBinError(i, 0)
        ratios = []
        for s in reversed(samples_sys):
            if (len(samples) > 1 or options.sherpa_pdf) and "_GEN_" in s:
                continue
            h_pt_truth_r = h_pt_truth[s].Clone(f"{h_pt_truth[s].GetName()}_r")
            h_pt_truth_r.Divide(denum)
            h_pt_truth_r.Draw("hist same")
            ratios += [h_pt_truth_r]

        # save
        ROOT.gPad.RedrawAxis()
        canv3.set_ratio_range(0.71, 1.29, override=True)
        canv3.print(f"{options.output}/{c}_truth.pdf")
        h_pt_truth[samples[0]].Write()

        # -------------------
        # draw fiducial efficiency
        # -------------------
        if options.sherpa_pdf or options.sherpa_qcd or options.sherpa_as or options.sherpa_ew:
            sherpa_nominal = None
            sherpa_sys = []
            for key, gr in fid_eff_gr.items():
                if "Sherpa2211_WplusD" not in key:
                    continue
                if key == "Sherpa2211_WplusD":
                    sherpa_nominal = gr
                else:
                    sherpa_sys += [gr]
            if options.sherpa_pdf:
                sys_band, sys_band_ratio = utils.make_pdf_err(sherpa_nominal, sherpa_sys, "NNPDF30_nnlo_as_0118_hessian")
            elif options.sherpa_qcd or options.sherpa_as or options.sherpa_ew:
                sys_band, sys_band_ratio = utils.make_minmax_err(sherpa_nominal, sherpa_sys)
        ROOT.gStyle.SetPaintTextFormat(".5f")

        # canvas
        canv4 = utils.make_canvas_mc_ratio(proxy_axis[samples[0]], truth_pt, chan, "Ratio", x=800, y=800, events="fiducial efficiency")
        canv4.pad1.cd()
        canv4.pad1.SetLogx()
        # if len(samples) < 2:
        #     for j, s in enumerate(samples):
        #         for i in range(1, nbins + 2):
        #             proxy_axis[s].SetBinContent(i, fid_eff_gr[s].GetY()[i - 1])
        #             proxy_axis[s].SetBinError(i, 0)
        #             proxy_axis[s].Draw("hist text same")
        # eff = fid_eff_inclusive[s].GetEfficiency(1)
        # eff_err = (fid_eff_inclusive[s].GetEfficiencyErrorUp(1) + fid_eff_inclusive[s].GetEfficiencyErrorLow(1)) / 2
        # ROOT.myText(0.18, 0.72 - 0.10 * (j + 1), 1, f"{s}: {eff:1.5f} #pm {eff_err:1.5f} ({100 * eff_err / eff:1.3f}%)")

        if options.sherpa_pdf or options.sherpa_qcd or options.sherpa_as or options.sherpa_ew:
            sys_band.SetFillColor(colors["Sherpa2211_WplusD"] + 2)
            sys_band.Draw("e2")
        for s in samples_sys:
            if (len(samples) > 1 or options.sherpa_pdf) and "_GEN_" in s:
                continue
            fid_eff_gr[s].Draw("pe")
        leg.Draw()

        # ratio
        canv4.pad2.cd()
        canv4.pad2.SetLogx()
        if len(samples) > 1:
            canv4.set_ratio_range(0.81, 1.19, override=True)
        else:
            canv4.set_ratio_range(0.951, 1.049, override=True)
        fid_eff_gr_ratio = {s: fid_eff_gr[s].Clone(f"{fid_eff_gr[s].GetName()}_ratio") for s in samples_sys}
        for s in samples_sys:
            for i in range(fid_eff_gr_ratio[s].GetN()):
                fid_eff_gr_ratio[s].GetY()[i] /= fid_eff_gr[samples[0]].GetY()[i]
                fid_eff_gr_ratio[s].GetEYhigh()[i] /= fid_eff_gr[samples[0]].GetY()[i]
                fid_eff_gr_ratio[s].GetEYlow()[i] /= fid_eff_gr[samples[0]].GetY()[i]

            if "_GEN_" in s:
                continue

            h, _, _ = utils.get_hist_from_gr(fid_eff_gr_ratio[s], f"{s}_{c}_fid_eff_ratio")
            h.Write()

        # draw Sherpa errors
        if options.sherpa_pdf or options.sherpa_qcd or options.sherpa_as or options.sherpa_ew:

            # adjust the central values
            sys_band_ratio_clone = sys_band_ratio.Clone()
            for i in range(fid_eff_gr_ratio["Sherpa2211_WplusD"].GetN()):
                sys_band_ratio_clone.GetY()[i] = fid_eff_gr_ratio["Sherpa2211_WplusD"].GetY()[i]
            sys_band_ratio_clone.Draw("e2")
            sys_band_ratio_clone.SetFillColor(colors["Sherpa2211_WplusD"] + 2)

            if options.sherpa_qcd or options.sherpa_as or options.sherpa_ew:
                _, h_up, h_dn = utils.get_hist_from_gr(sys_band_ratio, f"{s}_{c}_fid_eff_ratio_qcd_err")
                sys_band_ratio.Write(f"gr_{s}_{c}_fid_eff_ratio_qcd_err")
                h_up.Write()
                h_dn.Write()
            elif options.sherpa_pdf:
                _, h_up, h_dn = utils.get_hist_from_gr(sys_band_ratio, f"{s}_{c}_fid_eff_ratio_pdf_err")
                sys_band_ratio.Write(f"gr_{s}_{c}_fid_eff_ratio_pdf_err")
                h_up.Write()
                h_dn.Write()
        else:
            for s in samples_sys:
                h, _, _ = utils.get_hist_from_gr(fid_eff_gr_ratio[s], f"{s}_{c}_fid_eff_ratio")
                h.Write()

        # draw all graphs
        for s in reversed(samples_sys):
            if (len(samples) > 1 or options.sherpa_pdf) and "_GEN_" in s:
                continue
            fid_eff_gr_ratio[s].Draw("pe")

        # save
        ROOT.gPad.RedrawAxis()
        canv4.print(f"{options.output}/{c}_fid_eff.pdf")

        # -------------------
        # draw fiducial efficiency per bin
        # -------------------
        ROOT.gStyle.SetPaintTextFormat(".2f%")
        ROOT.gStyle.SetPalette(ROOT.kCMYK)
        canv5 = ROOT.TCanvas(f"{c}_matrix", f"{c}_matrix", 1000, 800)
        canv5.SetRightMargin(0.20)
        canv5.SetLogy()
        canv5.SetLogx()
        h_fid_eff[samples[0]].GetXaxis().SetMoreLogLabels()
        h_fid_eff[samples[0]].GetXaxis().SetNoExponent()
        h_fid_eff[samples[0]].GetYaxis().SetMoreLogLabels()
        h_fid_eff[samples[0]].GetZaxis().SetTitle("Response Matrix [%] #pm rel. err. [%]")
        h_fid_eff[samples[0]].Draw("text colz error")

        # ATLAS label
        ROOT.ATLASLabel(0.18, 0.90, "Internal", 1)
        ROOT.myText(0.18, 0.84, 1, "#sqrt{s} = 13 TeV")
        ROOT.myText(0.18, 0.78, 1, "139 fb^{-1}")
        if "Dplus" in c:
            ROOT.myText(0.40, 0.18, 1, "W#rightarrowl#nu+D, D#rightarrowK#pi#pi")
        elif "Dstar" in c:
            ROOT.myText(0.40, 0.18, 1, "W#rightarrowl#nu+D*, D*#rightarrow(K#pi)#pi")
        # ROOT.myText(0.50, 0.18, 1, c.replace("OS-SS_", ""))

        # save
        ROOT.gPad.RedrawAxis()
        canv5.Print(f"{options.output}/{c}_fid_eff_per_bin.pdf")
        h_fid_eff[samples[0]].Write()

        # normalized matrix
        integral = h_fid_eff[samples[0]].GetSumOfWeights()
        for i in range(1, h_fid_eff[samples[0]].GetNbinsX() + 1):
            for j in range(1, h_fid_eff[samples[0]].GetNbinsX() + 1):
                z = h_fid_eff[samples[0]].GetBinContent(i, j)
                h_fid_eff[samples[0]].SetBinContent(i, j, 100 * z / integral)
        h_fid_eff[samples[0]].SetMaximum(25.)
        h_fid_eff[samples[0]].GetZaxis().SetTitle("Normalized R. Matrix [%] #pm rel. err. [%]")
        canv5.Update()
        canv5.Print(f"{options.output}/{c}_fid_eff_per_bin_norm.pdf")

        # -------------------
        # draw the invertex matrix
        # -------------------
        ROOT.gStyle.SetPaintTextFormat(".3f")
        canv6 = ROOT.TCanvas(f"{c}_matrix_inv", f"{c}_matrix_inv", 1000, 800)
        canv6.SetRightMargin(0.15)
        canv6.SetLogy()
        canv6.SetLogx()
        h_fid_eff_inv[samples[0]].GetXaxis().SetMoreLogLabels()
        h_fid_eff_inv[samples[0]].GetXaxis().SetNoExponent()
        h_fid_eff_inv[samples[0]].GetYaxis().SetMoreLogLabels()
        h_fid_eff_inv[samples[0]].Draw("text colz")

        # ATLAS label
        # ROOT.ATLASLabel(0.18, 0.90, "Internal", 1)
        # ROOT.myText(0.18, 0.84, 1, "#sqrt{s} = 13 TeV")
        # ROOT.myText(0.18, 0.78, 1, "139 fb^{-1}")
        # ROOT.myText(0.50, 0.24, 1, "W#rightarrowl#nu+D, D#rightarrowK#pi#pi")
        # ROOT.myText(0.50, 0.18, 1, c.replace("OS-SS_", ""))

        # save
        ROOT.gPad.RedrawAxis()
        canv6.Print(f"{options.output}/{c}_fid_eff_inv.pdf")

    # close file
    f_out.Close()


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-o', '--output',
                      action="store", dest="output",
                      default="transfer_matrix", help="Output folder path.")
    parser.add_option('-t', '--truth',
                      action="store", dest="truth",
                      default="", help="Path to the histograms.root file from the truth analysis.")
    parser.add_option('-r', '--reco',
                      action="store", dest="reco",
                      default="", help="Path to the histograms.root file from the reco analysis.")
    parser.add_option('-s', '--samples',
                      action="store", dest="samples",
                      default="MG_Wjets", help="Samples to use. Comma separates list.")
    parser.add_option('-c', '--channels',
                      action="store", dest="channels",
                      default="OS-SS_Dplus_Kpipi", help="Channels to use. Comma separates list.")
    parser.add_option('--sherpa-pdf',
                      action="store_true", dest="sherpa_pdf",
                      default=False, help="Sherpa PDF systematics.")
    parser.add_option('--sherpa-qcd',
                      action="store_true", dest="sherpa_qcd",
                      default=False, help="Sherpa QCD systematics.")
    parser.add_option('--sherpa-as',
                      action="store_true", dest="sherpa_as",
                      default=False, help="Sherpa Alpha_S systematics.")
    parser.add_option('--sherpa-ew',
                      action="store_true", dest="sherpa_ew",
                      default=False, help="Sherpa EW systematics.")

    # parse input arguments
    options, args = parser.parse_args()

    # run
    main(options, args)
