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

# pT(D) bins
# last bin extends to infinity
# ends at 120 for cosmetical purposes
PT_BINS = [8, 12, 20, 40, 80, 120]

# max y for differential cross sections [pb]
Y_MAX = 75

# theory prediction style
THEORY_DICT = {
    "MG_Wjets": {
        "lineColor": ROOT.kBlue,
        "fillColor": ROOT.kBlue - 9,
        "markerStyle": 32,
        "legendLabel": "LO MG",
        "offset": -0.70,
    },
    "MGPy8EG_NLO_WplusD": {
        "lineColor": ROOT.kRed,
        "fillColor": ROOT.kRed - 9,
        "markerStyle": 26,
        "legendLabel": "NLO MG #it{W}+#it{D}",
        "offset": -0.35,
    },
    "Sherpa2211_WplusD": {
        "lineColor": ROOT.kCyan + 2,
        "fillColor": ROOT.kCyan - 4,
        "markerStyle": 27,
        "legendLabel": "Sh2.2.11 #it{W}+#it{D}",
        "offset": 0.35,
    },
    "MGFxFx_WplusD": {
        "lineColor": ROOT.kMagenta + 1,
        "fillColor": ROOT.kMagenta - 9,
        "markerStyle": 4,
        "legendLabel": "MG FxFx #it{W}+#it{D}",
        "offset": 0.70,
    },
}

# make output folder
if not os.path.isdir("fit_results"):
    os.makedirs("fit_results")


# extract POIs from txt file
def extract_pois(file):
    POIs = {}
    with open(file) as f:
        for line in f:
            res = line.split()
            if not len(res):
                continue
            name = res[0]
            if "mu_" in name:
                POIs[name] = res[1:]
    return POIs


def createCanvasPads(name):
    c = ROOT.TCanvas(name, name, 1200, 1200)
    # Upper histogram plot is pad1
    pad1 = ROOT.TPad(f"pad1_{name}", f"pad1_{name}", 0, 0.45, 1, 1.0)
    pad1.SetTopMargin(0.1)
    pad1.SetBottomMargin(0.05 / 0.55)
    pad1.SetFillStyle(4000)
    pad1.Draw()
    # Middle ratio plot is pad2
    c.cd()
    pad2 = ROOT.TPad(f"pad2_{name}", f"pad2_{name}", 0, 0.20, 1, 0.55)
    pad2.SetTopMargin(0.05 / 0.35)
    pad2.SetBottomMargin(0.05 / 0.35)
    pad2.SetFillStyle(4000)
    pad2.Draw()
    # Lower ratio plot is pad3
    c.cd()
    pad3 = ROOT.TPad(f"pad3_{name}", f"pad3_{name}", 0, 0.00, 1, 0.30)
    pad3.SetTopMargin(0.05 / 0.30)
    pad3.SetBottomMargin(0.10 / 0.30)
    pad3.SetFillStyle(4000)
    pad3.Draw()

    return c, pad1, pad2, pad3


# --------------------------------------------
# Step 0: get cross section priors
# --------------------------------------------
dir_prior = "/global/cfs/cdirs/atlas/wcharm/charmplot_output/Dmeson_2022_03_15/"
f = ROOT.TFile(os.path.join(dir_prior, "fid_eff_dplus_stat", "unfolding.root"))
h_minus = f.Get("Sherpa2211_WplusD_OS-SS_lep_plus_0tag_Dplus_Kpipi_truth_differential_pt")
h_plus = f.Get("Sherpa2211_WplusD_OS-SS_lep_plus_0tag_Dplus_Kpipi_truth_differential_pt")
priors = {
    "Wminus": h_minus.Integral(),
    "Wminus_1": h_minus.GetBinContent(1),
    "Wminus_2": h_minus.GetBinContent(2),
    "Wminus_3": h_minus.GetBinContent(3),
    "Wminus_4": h_minus.GetBinContent(4),
    "Wminus_5": h_minus.GetBinContent(5),
    "Wplus": h_plus.Integral(),
    "Wplus_1": h_plus.GetBinContent(1),
    "Wplus_2": h_plus.GetBinContent(2),
    "Wplus_3": h_plus.GetBinContent(3),
    "Wplus_4": h_plus.GetBinContent(4),
    "Wplus_5": h_plus.GetBinContent(5),
}

# --------------------------------------------
# Step 1: parse fit results from txt files
# --------------------------------------------
# folder names
stat_only = "WCharm_lep_obs_stat_only_OSSS_complete"
obs = "WCharm_lep_obs_OSSS_complete"

# input path
dir = "/global/cfs/cdirs/atlas/wcharm/TRExFitter/Output/Dplus_2022_03_07/WCharm_lep_OSSS_complete_all"

# stat-only
POIs_stat = extract_pois(os.path.join(dir, stat_only, "Fits", f"{stat_only}.txt"))
POIs_stat.update(extract_pois(os.path.join(dir, stat_only, "Fits", f"{stat_only}_expr.txt")))

# observed
POIs_obs = extract_pois(os.path.join(dir, obs, "Fits", f"{obs}.txt"))
POIs_obs.update(extract_pois(os.path.join(dir, obs, "Fits", f"{obs}_expr.txt")))
print(POIs_obs)

# --------------------------------------------
# Step 2: make TGraph objects
# --------------------------------------------
trash = []
for lep in ["minus", "plus"]:
    gr_obs = ROOT.TGraphAsymmErrors()
    gr_obs_norm = ROOT.TGraphAsymmErrors()
    gr_obs_ratio = ROOT.TGraphAsymmErrors()
    gr_obs_norm_ratio = ROOT.TGraphAsymmErrors()
    gr_obs_norm_ratio_stat = ROOT.TGraphAsymmErrors()
    gr_obs_sys = ROOT.TGraphAsymmErrors()
    gr_obs_norm_sys = ROOT.TGraphAsymmErrors()
    for i in range(len(PT_BINS) - 1):
        xl = PT_BINS[i]
        xh = PT_BINS[i + 1]
        w = xh - xl

        y_prior = float(priors[f"W{lep}_{i + 1}"])
        y_rel = y_prior / float(priors[f"W{lep}"])
        expr = ""
        if lep == "plus":
            expr = "expr_"
        y_tot = float(POIs_obs[f"{expr}mu_W{lep}_tot"][0]) * y_prior

        if i < 4:
            y = float(POIs_obs[f"expr_mu_W{lep}_{i + 1}"][0]) * y_prior
            y_up = float(POIs_obs[f"expr_mu_W{lep}_{i + 1}"][1]) * y_prior
            y_dn = float(POIs_obs[f"expr_mu_W{lep}_{i + 1}"][1]) * y_prior
            y_up_stat = float(POIs_stat[f"expr_mu_W{lep}_{i + 1}"][1]) * y_prior
            y_dn_stat = float(POIs_stat[f"expr_mu_W{lep}_{i + 1}"][1]) * y_prior
            y_up_sys = (y_up**2 - y_up_stat**2)**(0.5)
            y_dn_sys = (y_dn**2 - y_dn_stat**2)**(0.5)
            y_norm = float(POIs_obs[f"mu_W{lep}_rel_{i + 1}"][0]) * y_rel
            y_norm_up = float(POIs_obs[f"mu_W{lep}_rel_{i + 1}"][1]) * y_rel
            y_norm_dn = float(POIs_obs[f"mu_W{lep}_rel_{i + 1}"][2]) * y_rel
            y_norm_up_stat = float(POIs_stat[f"mu_W{lep}_rel_{i + 1}"][1]) * y_rel
            y_norm_dn_stat = float(POIs_stat[f"mu_W{lep}_rel_{i + 1}"][2]) * y_rel
        else:
            y = float(POIs_obs[f"expr_mu_W{lep}_rel_{i + 1}"][0]) * y_tot
            y_up = ((float(POIs_obs[f"expr_mu_W{lep}_rel_{i + 1}"][1]) * y_tot)**2 + (float(POIs_obs[f"{expr}mu_W{lep}_tot"][1]) * y_prior)**2)**(0.5)
            y_dn = ((float(POIs_obs[f"expr_mu_W{lep}_rel_{i + 1}"][1]) * y_tot)**2 + (float(POIs_obs[f"{expr}mu_W{lep}_tot"][1]) * y_prior)**2)**(0.5)
            y_up_stat = ((float(POIs_stat[f"expr_mu_W{lep}_rel_{i + 1}"][1]) * y_tot)**2 + (float(POIs_stat[f"{expr}mu_W{lep}_tot"][1]) * y_prior)**2)**(0.5)
            y_dn_stat = ((float(POIs_stat[f"expr_mu_W{lep}_rel_{i + 1}"][1]) * y_tot)**2 + (float(POIs_stat[f"{expr}mu_W{lep}_tot"][1]) * y_prior)**2)**(0.5)
            y_up_sys = (y_up**2 - y_up_stat**2)**(0.5)
            y_dn_sys = (y_dn**2 - y_dn_stat**2)**(0.5)
            y_norm = float(POIs_obs[f"expr_mu_W{lep}_rel_{i + 1}"][0]) * y_rel
            y_norm_up = float(POIs_obs[f"expr_mu_W{lep}_rel_{i + 1}"][1]) * y_rel
            y_norm_dn = float(POIs_obs[f"expr_mu_W{lep}_rel_{i + 1}"][1]) * y_rel
            y_norm_up_stat = float(POIs_stat[f"expr_mu_W{lep}_rel_{i + 1}"][1]) * y_rel
            y_norm_dn_stat = float(POIs_stat[f"expr_mu_W{lep}_rel_{i + 1}"][1]) * y_rel
        y_norm_up_sys = (y_norm_up**2 - y_norm_up_stat**2)**(0.5)
        y_norm_dn_sys = (y_norm_dn**2 - y_norm_dn_stat**2)**(0.5)

        # fill graphs
        # xc = xl + w / 2.
        xc = ROOT.TMath.Power(10, ROOT.TMath.Log10(xl) + (ROOT.TMath.Log10(xl + w) - ROOT.TMath.Log10(xl)) / 2.)
        xc_up = xl + w - xc
        xc_dn = xc - xl
        gr_obs.SetPoint(i, xc, y)
        gr_obs.SetPointError(i, xc_dn, xc_up, abs(y_dn), abs(y_up))
        gr_obs_norm.SetPoint(i, xc, y_norm)
        gr_obs_norm.SetPointError(i, xc_dn, xc_up, abs(y_norm_dn), abs(y_norm_up))
        gr_obs_ratio.SetPoint(i, xc, 1.0)
        gr_obs_ratio.SetPointError(i, xc_dn, xc_up, abs(y_dn / y), abs(y_up / y))
        gr_obs_norm_ratio.SetPoint(i, xc, 1.0)
        gr_obs_norm_ratio.SetPointError(i, xc_dn, xc_up, abs(y_norm_dn / y_norm), abs(y_norm_up / y_norm))
        gr_obs_norm_ratio_stat.SetPoint(i, xc, 1.0)
        gr_obs_norm_ratio_stat.SetPointError(i, xc_dn, xc_up, abs(y_norm_dn_stat / y_norm), abs(y_norm_up_stat / y_norm))
        gr_obs_sys.SetPoint(i, xc, y)
        gr_obs_sys.SetPointError(i, w / 10., w / 10., abs(y_dn_sys), abs(y_up_sys))
        gr_obs_norm_sys.SetPoint(i, xc, y_norm)
        gr_obs_norm_sys.SetPointError(i, w / 10., w / 10., abs(y_norm_dn_sys), abs(y_norm_up_sys))

    # line width
    gr_obs.SetLineWidth(2)
    gr_obs_norm.SetLineWidth(2)
    gr_obs_ratio.SetLineWidth(0)
    gr_obs_norm_ratio.SetLineWidth(0)
    gr_obs.SetMarkerSize(1.5)
    gr_obs_norm.SetMarkerSize(1.5)

    # fill
    gr_obs_sys.SetLineColor(ROOT.kBlack)
    gr_obs_norm_sys.SetLineColor(ROOT.kBlack)
    gr_obs_ratio.SetFillColor(ROOT.kGray)
    # gr_obs_norm_ratio.SetFillColor(ROOT.kGray + 1)
    gr_obs_norm_ratio.SetFillColor(ROOT.kGray)
    gr_obs_norm_ratio_stat.SetFillColor(ROOT.kGray + 1)
    gr_obs_norm_ratio_stat.SetLineColor(ROOT.kGray + 1)

    # multigraph
    mg_obs = ROOT.TMultiGraph()
    mg_obs_norm = ROOT.TMultiGraph()
    mg_obs_ratio = ROOT.TMultiGraph()

    # axis title
    mg_obs.GetYaxis().SetTitle("d#sigma/d#it{p}_{T}^{#it{D}} [pb]")
    mg_obs_norm.GetYaxis().SetTitle("1/#sigma^{tot.}_{fid.} d#sigma/d#it{p}_{T}^{#it{D}}")
    mg_obs_ratio.GetYaxis().SetTitle("#frac{Theory}{1/#sigma^{tot.}_{fid.} d#sigma/d#it{p}_{T}^{#it{D}}}")
    mg_obs_ratio.GetXaxis().SetTitle("#it{p}_{T}(#it{D}) [GeV]")

    # label and title size
    GLOBAL_SF = 1.3
    mg_obs.GetYaxis().SetTitleSize(mg_obs.GetYaxis().GetTitleSize() * GLOBAL_SF)
    mg_obs.GetYaxis().SetTitleOffset(mg_obs.GetYaxis().GetTitleOffset() * (1 / (GLOBAL_SF * 1.1)))
    mg_obs.GetYaxis().SetLabelSize(mg_obs.GetYaxis().GetLabelSize() * GLOBAL_SF)

    SF = 0.55 / 0.35
    mg_obs_norm.GetYaxis().SetTitleSize(mg_obs_norm.GetYaxis().GetTitleSize() * SF * GLOBAL_SF)
    mg_obs_norm.GetYaxis().SetTitleOffset(mg_obs_norm.GetYaxis().GetTitleOffset() * (1 / (GLOBAL_SF * 1.1 * SF)))
    mg_obs_norm.GetYaxis().SetLabelSize(mg_obs_norm.GetYaxis().GetLabelSize() * SF * GLOBAL_SF)

    SF = 0.55 / 0.30
    mg_obs_ratio.GetYaxis().SetTitleSize(mg_obs_ratio.GetYaxis().GetTitleSize() * SF * GLOBAL_SF)
    mg_obs_ratio.GetYaxis().SetTitleOffset(mg_obs_ratio.GetYaxis().GetTitleOffset() * (1 / (GLOBAL_SF * 1.0 * SF)))
    mg_obs_ratio.GetXaxis().SetTitleSize(mg_obs_ratio.GetXaxis().GetTitleSize() * SF * GLOBAL_SF)
    mg_obs_ratio.GetXaxis().SetTitleOffset(mg_obs_ratio.GetXaxis().GetTitleOffset() * (1 / (GLOBAL_SF * 0.5 * SF)))
    mg_obs_ratio.GetYaxis().SetLabelSize(mg_obs_ratio.GetYaxis().GetLabelSize() * SF * GLOBAL_SF)
    mg_obs_ratio.GetXaxis().SetLabelSize(mg_obs_ratio.GetXaxis().GetLabelSize() * SF * GLOBAL_SF)

    # tick marks
    mg_obs.GetYaxis().SetNdivisions(506)
    mg_obs_norm.GetYaxis().SetNdivisions(504)
    mg_obs_ratio.GetYaxis().SetNdivisions(306)

    # log x-axis
    mg_obs_ratio.GetXaxis().SetMoreLogLabels()
    mg_obs_ratio.GetXaxis().SetNoExponent()

    # legend
    N = 3 + len(THEORY_DICT)
    leg = ROOT.TLegend(0.66, 0.87 - N * (0.067), 0.98, 0.87)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(36)
    leg.SetTextFont(43)
    leg.AddEntry(gr_obs, "Data", "pe")
    leg.AddEntry(gr_obs_sys, "Syst. Unc.", "f")
    leg.AddEntry(gr_obs_norm_ratio_stat, "Stat. Unc.", "f")
    leg.AddEntry(gr_obs_norm_ratio, "Syst. #oplus Stat.", "f")

    # add to ratio multigraph
    mg_obs_ratio.Add(gr_obs_norm_ratio, "e2")
    mg_obs_ratio.Add(gr_obs_norm_ratio_stat, "e2")

    # --------------------------------------------
    # Step 2.5: theory comparisons
    # --------------------------------------------
    for prediction in ["MG_Wjets", "MGPy8EG_NLO_WplusD", "Sherpa2211_WplusD", "MGFxFx_WplusD"]:
        f_theory = ROOT.TFile(os.path.join(dir_prior, "histograms.root"))
        h = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_D_pt_fit")
        gr_qcd = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_D_pt_fit_ratio_sherpa2211_theory_qcd")
        gr_pdf = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_D_pt_fit_ratio_sherpa2211_theory_pdf")
        gr_as = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_D_pt_fit_ratio_sherpa2211_theory_as")
        gr_theory = gr_obs.Clone()
        total_xsec = 0.
        total_xsec_up = 0.
        total_xsec_dn = 0.
        BR = 1.0
        # if prediction == "MG_Wjets":
        #     BR = 0.094
        for i in range(h.GetNbinsX()):
            gr_theory.GetY()[i] = h.GetBinContent(i + 1) / BR
            if THEORY_DICT[prediction]["offset"] > 0:
                gr_theory.GetX()[i] = gr_theory.GetX()[i] + THEORY_DICT[prediction]["offset"] * gr_theory.GetEXhigh()[i]
            else:
                gr_theory.GetX()[i] = gr_theory.GetX()[i] + THEORY_DICT[prediction]["offset"] * gr_theory.GetEXlow()[i]
            gr_theory.GetEXhigh()[i] = gr_theory.GetEXhigh()[i] / 10.
            gr_theory.GetEXlow()[i] = gr_theory.GetEXlow()[i] / 10.
            gr_theory.GetEYhigh()[i] = ((h.GetBinError(i + 1) / BR)**2 +
                                        (gr_qcd.GetEYhigh()[i] * gr_theory.GetY()[i])**2 +
                                        (gr_pdf.GetEYhigh()[i] * gr_theory.GetY()[i])**2 +
                                        (gr_as.GetEYhigh()[i] * gr_theory.GetY()[i])**2)**(0.5)
            gr_theory.GetEYlow()[i] = ((h.GetBinError(i + 1) / BR)**2 +
                                       (gr_qcd.GetEYlow()[i] * gr_theory.GetY()[i])**2 +
                                       (gr_pdf.GetEYlow()[i] * gr_theory.GetY()[i])**2 +
                                       (gr_as.GetEYlow()[i] * gr_theory.GetY()[i])**2)**(0.5)
            total_xsec += gr_theory.GetY()[i]
            total_xsec_up += gr_theory.GetY()[i] + gr_theory.GetEYhigh()[i]
            total_xsec_dn += gr_theory.GetY()[i] - gr_theory.GetEYlow()[i]

        # style
        gr_theory.SetLineColor(THEORY_DICT[prediction]["lineColor"])
        gr_theory.SetMarkerColor(THEORY_DICT[prediction]["lineColor"])
        gr_theory.SetFillColor(THEORY_DICT[prediction]["fillColor"])
        gr_theory.SetMarkerStyle(THEORY_DICT[prediction]["markerStyle"])

        # normalized cross section
        gr_theory_norm = gr_theory.Clone()
        for i in range(h.GetNbinsX()):
            gr_theory_norm.GetY()[i] = gr_theory_norm.GetY()[i] / total_xsec
            gr_theory_norm.GetEYhigh()[i] = (gr_theory.GetY()[i] + gr_theory.GetEYhigh()[i]) / total_xsec_up - gr_theory_norm.GetY()[i]
            gr_theory_norm.GetEYlow()[i] = gr_theory_norm.GetY()[i] - (gr_theory.GetY()[i] - gr_theory.GetEYlow()[i]) / total_xsec_dn

        # ratio plot
        gr_theory_ratio = gr_theory_norm.Clone()
        for i in range(h.GetNbinsX()):
            gr_theory_ratio.GetY()[i] = gr_theory_norm.GetY()[i] / gr_obs_norm.GetY()[i]
            gr_theory_ratio.GetEYhigh()[i] = gr_theory_norm.GetEYhigh()[i] / gr_obs_norm.GetY()[i]
            gr_theory_ratio.GetEYlow()[i] = gr_theory_norm.GetEYlow()[i] / gr_obs_norm.GetY()[i]

        # add to legend
        leg.AddEntry(gr_theory_ratio, THEORY_DICT[prediction]["legendLabel"], "pf")

        # add to multigraph
        mg_obs.Add(gr_theory, "pe5")
        mg_obs_norm.Add(gr_theory_norm, "pe5")
        mg_obs_ratio.Add(gr_theory_ratio, "pe5")

    # --------------------------------------------
    # Step 3: make plots
    # --------------------------------------------
    # add graphs
    mg_obs.Add(gr_obs_sys, "e5")
    mg_obs.Add(gr_obs, "pe")
    mg_obs_norm.Add(gr_obs_norm_sys, "e5")
    mg_obs_norm.Add(gr_obs_norm, "pe")

    # Plot histograms: first canvas, draw OS and SS on the same plot
    c1, pad1, pad2, pad3 = createCanvasPads(f"W{lep}")

    # upper canvas
    pad1.cd()
    pad1.SetLogx()
    mg_obs.Draw("a")
    mg_obs.GetXaxis().SetLimits(8, 120)
    mg_obs.SetMinimum(1e-3)
    mg_obs.SetMaximum(Y_MAX)

    # ATLAS label
    l1 = ROOT.TLatex()
    l1.SetTextFont(73)
    l1.SetTextSize(36)
    l1.DrawLatex(9, Y_MAX * (1 - 0.10), "ATLAS")
    l2 = ROOT.TLatex()
    l2.SetTextFont(43)
    l2.SetTextSize(36)
    l2.DrawLatex(13, Y_MAX * (0.9 - 0 * 0.08), "Internal")
    l2.DrawLatex(9, Y_MAX * (0.9 - 1 * 0.08), "#sqrt{s} = 13 TeV, 139 fb^{-1}")
    l2.DrawLatex(9, Y_MAX * (0.9 - 2 * 0.08), "#it{W}+#it{D}(#rightarrowK#pi#pi), #it{W%s} channel" % ("-" if lep == "minus" else "+"))

    # legend
    leg.Draw()

    ROOT.gPad.RedrawAxis()
    pad2.cd()
    pad2.SetLogx()
    # pad2.SetLogy()
    mg_obs_norm.Draw("a")
    mg_obs_norm.GetXaxis().SetLimits(8, 120)
    # mg_obs_norm.SetMinimum(0.015)
    mg_obs_norm.SetMinimum(1e-3)
    mg_obs_norm.SetMaximum(0.35)

    ROOT.gPad.RedrawAxis()
    pad3.cd()
    pad3.SetLogx()
    mg_obs_ratio.Draw("a")
    mg_obs_ratio.GetXaxis().SetLimits(8, 120)
    mg_obs_ratio.SetMinimum(0.83)
    mg_obs_ratio.SetMaximum(1.17)
    line3 = ROOT.TLine(8, 1, 120, 1)
    line3.SetLineStyle(2)
    line3.Draw()

    ROOT.gPad.RedrawAxis()
    c1.Print(f"fit_results/W{lep}.pdf")
    c1.Print(f"fit_results/W{lep}.png")

    trash += [mg_obs]
    trash += [mg_obs_norm]
    trash += [mg_obs_ratio]

f.Close()
