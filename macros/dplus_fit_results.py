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
Y_MAX = 50

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
dir = "/global/cfs/cdirs/atlas/wcharm/charmplot_output/Dmeson_2022_02_02/fid_eff_dplus_stat"
f = ROOT.TFile(os.path.join(dir, "unfolding.root"))
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
            y_norm_up_sys = (y_norm_up**2 - (float(POIs_stat[f"mu_W{lep}_rel_{i + 1}"][1]) * y_rel)**2)**(0.5)
            y_norm_dn_sys = (y_norm_dn**2 - (float(POIs_stat[f"mu_W{lep}_rel_{i + 1}"][1]) * y_rel)**2)**(0.5)
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
            y_norm_up_sys = (y_norm_up**2 - (float(POIs_stat[f"expr_mu_W{lep}_rel_{i + 1}"][1]) * y_rel)**2)**(0.5)
            y_norm_dn_sys = (y_norm_dn**2 - (float(POIs_stat[f"expr_mu_W{lep}_rel_{i + 1}"][1]) * y_rel)**2)**(0.5)

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
    gr_obs_norm_ratio.SetFillColor(ROOT.kGray + 1)

    # multigraph
    mg_obs = ROOT.TMultiGraph()
    mg_obs_norm = ROOT.TMultiGraph()
    mg_obs_ratio = ROOT.TMultiGraph()

    # add graphs
    mg_obs.Add(gr_obs_sys, "e5")
    mg_obs.Add(gr_obs, "pe")
    mg_obs_norm.Add(gr_obs_norm_sys, "e5")
    mg_obs_norm.Add(gr_obs_norm, "pe")
    mg_obs_ratio.Add(gr_obs_ratio, "e2")
    mg_obs_ratio.Add(gr_obs_norm_ratio, "e2")

    # axis title
    mg_obs.GetYaxis().SetTitle("d#sigma/d#it{p}_{T}^{#it{D}} [pb]")
    mg_obs_norm.GetYaxis().SetTitle("1/#sigma^{tot.}_{fid.} d#sigma/d#it{p}_{T}^{#it{D}}")
    mg_obs_ratio.GetYaxis().SetTitle("Theory / Data")
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
    mg_obs_ratio.GetYaxis().SetTitleOffset(mg_obs_ratio.GetYaxis().GetTitleOffset() * (1 / (GLOBAL_SF * 1.1 * SF)))
    mg_obs_ratio.GetXaxis().SetTitleSize(mg_obs_ratio.GetXaxis().GetTitleSize() * SF * GLOBAL_SF)
    mg_obs_ratio.GetXaxis().SetTitleOffset(mg_obs_ratio.GetXaxis().GetTitleOffset() * (1 / (GLOBAL_SF * 0.5 * SF)))
    mg_obs_ratio.GetYaxis().SetLabelSize(mg_obs_ratio.GetYaxis().GetLabelSize() * SF * GLOBAL_SF)
    mg_obs_ratio.GetXaxis().SetLabelSize(mg_obs_ratio.GetXaxis().GetLabelSize() * SF * GLOBAL_SF)

    # tick marks
    mg_obs.GetYaxis().SetNdivisions(506)
    mg_obs_norm.GetYaxis().SetNdivisions(503)
    mg_obs_ratio.GetYaxis().SetNdivisions(306)

    # log x-axis
    mg_obs_ratio.GetXaxis().SetMoreLogLabels()
    mg_obs_ratio.GetXaxis().SetNoExponent()

    # legend
    N = 4
    leg = ROOT.TLegend(0.60, 0.87 - N * (0.067), 0.95, 0.87)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(36)
    leg.SetTextFont(43)
    leg.AddEntry(gr_obs, "Data", "pe")
    leg.AddEntry(gr_obs_sys, "Syst. Unc.", "f")
    leg.AddEntry(gr_obs_ratio, "Syst. #oplus Stat. Unc.", "f")
    leg.AddEntry(gr_obs_norm_ratio, "1/#sigma^{tot.}_{fid.} Unc.", "f")

    # --------------------------------------------
    # Step 3: make plots
    # --------------------------------------------
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
    pad2.SetLogy()
    mg_obs_norm.Draw("a")
    mg_obs_norm.GetXaxis().SetLimits(8, 120)
    mg_obs_norm.SetMinimum(0.015)
    mg_obs_norm.SetMaximum(0.4)

    ROOT.gPad.RedrawAxis()
    pad3.cd()
    pad3.SetLogx()
    mg_obs_ratio.Draw("a")
    mg_obs_ratio.GetXaxis().SetLimits(8, 120)
    mg_obs_ratio.SetMinimum(0.91)
    mg_obs_ratio.SetMaximum(1.09)
    line3 = ROOT.TLine(8, 1, 120, 1)
    line3.SetLineStyle(2)
    line3.Draw()

    ROOT.gPad.RedrawAxis()
    c1.Print(f"fit_results/W{lep}.pdf")

    trash += [mg_obs]
    trash += [mg_obs_norm]
    trash += [mg_obs_ratio]
