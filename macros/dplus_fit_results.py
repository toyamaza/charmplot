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
trash = []

# cross section priors
DIR_PRIORS = "/global/cfs/cdirs/atlas/wcharm/charmplot_output/Dmeson_2022_05_07/"

# theory predictions
THEORY = ["MG_Wjets", "MGPy8EG_NLO_WplusD", "Sherpa2211_WplusD", "MGFxFx_WplusD"]

# max y for differential cross sections [pb]
Y_MAX = 75

# observables
OBSERVABLES = {
    "pt": {
        "fit_results": "/global/cfs/cdirs/atlas/wcharm/TRExFitter/Output/Dplus_2022_05_09_pt/",
        "label": "#it{p}_{T}^{#it{D}}",
        "prior_var": "D_pt_fit",
        "bins": [8, 12, 20, 40, 80, 120],
        "logx": True,
    },
    "eta": {
        "fit_results": "/global/cfs/cdirs/atlas/wcharm/TRExFitter/Output/Dplus_2022_05_09_eta/",
        "label": "|#eta(#it{l})|",
        "prior_var": "D_differential_lep_eta",
        "bins": [0.0, 0.5, 1.0, 1.5, 2.0, 2.5],
        "logx": False,
    },
}

# theory prediction style
THEORY_DICT = {
    "MG_Wjets": {
        "lineColor": ROOT.kBlue,
        "fillColor": ROOT.kBlue - 9,
        "markerStyle": 32,
        "markerStyleFull": 23,
        "legendLabel": "LO MG",
        "legendLabelFull": "LO MG",
        "offset": -0.70,
    },
    "MGPy8EG_NLO_WplusD": {
        "lineColor": ROOT.kRed,
        "fillColor": ROOT.kRed - 9,
        "markerStyle": 26,
        "markerStyleFull": 22,
        "legendLabel": "NLO MG #it{W}+#it{D}",
        "legendLabelFull": "NLO MG (NNPDF30_nnlo)",
        "offset": -0.35,
    },
    "Sherpa2211_WplusD": {
        "lineColor": ROOT.kGreen + 2,
        "fillColor": ROOT.kGreen - 6,
        "markerStyle": 27,
        "markerStyleFull": 33,
        "legendLabel": "Sh2.2.11 #it{W}+#it{D}",
        "legendLabelFull": "Sh2.2.11 (NNPDF30_nnlo)",
        "offset": 0.35,
    },
    "MGFxFx_WplusD": {
        "lineColor": ROOT.kMagenta + 1,
        "fillColor": ROOT.kMagenta - 9,
        "markerStyle": 4,
        "markerStyleFull": 20,
        "legendLabel": "MG FxFx #it{W}+#it{D}",
        "legendLabelFull": "MG FxFx (NNPDF31_nnlo)",
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


for obs_name, obs in OBSERVABLES.items():

    # no ticks in error bars
    ROOT.gStyle.SetEndErrorSize(0)

    # --------------------------------------------
    # Step 0: get cross section priors
    # --------------------------------------------
    f = ROOT.TFile(os.path.join(DIR_PRIORS, f"fid_eff_{obs_name}_dplus_stat", "unfolding.root"))
    h_minus = f.Get(f"Sherpa2211_WplusD_OS-SS_lep_minus_0tag_Dplus_Kpipi_truth_differential_{obs_name}")
    h_plus = f.Get(f"Sherpa2211_WplusD_OS-SS_lep_plus_0tag_Dplus_Kpipi_truth_differential_{obs_name}")
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
    print("============ cross section priors ============")
    for key, val in priors.items():
        print(f"{key}: {val}")

    # --------------------------------------------
    # Step 1: parse fit results from txt files
    # --------------------------------------------
    # folder names
    stat_only = "WCharm_lep_obs_stat_only_OSSS_complete"
    obs_fit = "WCharm_lep_obs_OSSS_complete"

    # stat-only
    POIs_stat = extract_pois(os.path.join(obs["fit_results"], stat_only, "Fits", f"{stat_only}.txt"))
    POIs_stat.update(extract_pois(os.path.join(obs["fit_results"], stat_only, "Fits", f"{stat_only}_expr.txt")))

    # observed
    POIs_obs = extract_pois(os.path.join(obs["fit_results"], obs_fit, "Fits", f"{obs_fit}.txt"))
    POIs_obs.update(extract_pois(os.path.join(obs["fit_results"], obs_fit, "Fits", f"{obs_fit}_expr.txt")))
    print("============ post fit results ============")
    for key, val in POIs_obs.items():
        print(f"{key}: {val}")

    # --------------------------------------------
    # Step 2: make TGraph objects
    # --------------------------------------------
    for lep in ["minus", "plus"]:
        gr_obs = ROOT.TGraphAsymmErrors()
        gr_obs_norm = ROOT.TGraphAsymmErrors()
        gr_obs_ratio = ROOT.TGraphAsymmErrors()
        gr_obs_norm_ratio = ROOT.TGraphAsymmErrors()
        gr_obs_norm_ratio_stat = ROOT.TGraphAsymmErrors()
        gr_obs_sys = ROOT.TGraphAsymmErrors()
        gr_obs_norm_sys = ROOT.TGraphAsymmErrors()
        w_tot = obs["bins"][-1] - obs["bins"][0]
        for i in range(len(obs["bins"]) - 1):
            xl = obs["bins"][i]
            xh = obs["bins"][i + 1]
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
                y_up_stat = ((float(POIs_stat[f"expr_mu_W{lep}_rel_{i + 1}"][1]) * y_tot)**2 +
                             (float(POIs_stat[f"{expr}mu_W{lep}_tot"][1]) * y_prior)**2)**(0.5)
                y_dn_stat = ((float(POIs_stat[f"expr_mu_W{lep}_rel_{i + 1}"][1]) * y_tot)**2 +
                             (float(POIs_stat[f"{expr}mu_W{lep}_tot"][1]) * y_prior)**2)**(0.5)
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
            if obs["logx"]:
                xc = ROOT.TMath.Power(10, ROOT.TMath.Log10(xl) + (ROOT.TMath.Log10(xl + w) - ROOT.TMath.Log10(xl)) / 2.)
            else:
                xc = xl + w / 2.
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
            gr_obs_norm_sys.SetPoint(i, xc, y_norm)
            if obs["logx"]:
                x_err_up = ROOT.TMath.Power(10, ROOT.TMath.Log10(xc) + (ROOT.TMath.Log10(obs["bins"][-1]) - ROOT.TMath.Log10(obs["bins"][0])) / 100.) - xc
                x_err_dn = xc - ROOT.TMath.Power(10, ROOT.TMath.Log10(xc) - (ROOT.TMath.Log10(obs["bins"][-1]) - ROOT.TMath.Log10(obs["bins"][0])) / 100.)
                gr_obs_sys.SetPointError(i, x_err_up, x_err_up, abs(y_dn_sys), abs(y_up_sys))
                gr_obs_norm_sys.SetPointError(i, x_err_up, x_err_up, abs(y_norm_dn_sys), abs(y_norm_up_sys))
            else:
                gr_obs_sys.SetPointError(i, xc_dn / 10., xc_up / 10., abs(y_dn_sys), abs(y_up_sys))
                gr_obs_norm_sys.SetPointError(i, xc_dn / 10., xc_up / 10., abs(y_norm_dn_sys), abs(y_norm_up_sys))

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
        gr_obs_ratio.SetFillColor(ROOT.kGray + 1)
        # gr_obs_norm_ratio.SetFillColor(ROOT.kGray + 1)
        gr_obs_norm_ratio.SetFillColor(ROOT.kGray + 1)
        gr_obs_norm_ratio_stat.SetFillColor(ROOT.kGray + 2)
        gr_obs_norm_ratio_stat.SetLineColor(ROOT.kGray + 2)

        # multigraph
        mg_obs = ROOT.TMultiGraph()
        mg_obs_norm = ROOT.TMultiGraph()
        mg_obs_ratio = ROOT.TMultiGraph()

        # axis title
        mg_obs.GetYaxis().SetTitle(f"d#sigma/d{obs['label']} [pb]")
        mg_obs_norm.GetYaxis().SetTitle(f"1/#sigma^{{tot.}}_{{fid.}} d#sigma/d{obs['label']}")
        mg_obs_ratio.GetYaxis().SetTitle(f"#frac{{Theory}}{{1/#sigma^{{tot.}}_{{fid.}} d#sigma/{obs['label']}}}")
        mg_obs_ratio.GetXaxis().SetTitle(obs['label'])

        # label and title size
        GLOBAL_SF = 1.3
        mg_obs.GetXaxis().SetLabelSize(0)
        mg_obs.GetYaxis().SetTitleSize(mg_obs.GetYaxis().GetTitleSize() * GLOBAL_SF)
        mg_obs.GetYaxis().SetTitleOffset(mg_obs.GetYaxis().GetTitleOffset() * (1 / (GLOBAL_SF * 1.1)))
        mg_obs.GetYaxis().SetLabelSize(mg_obs.GetYaxis().GetLabelSize() * GLOBAL_SF)

        SF = 0.55 / 0.35
        mg_obs_norm.GetXaxis().SetLabelSize(0)
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
        for prediction in THEORY:
            f_theory = ROOT.TFile(os.path.join(DIR_PRIORS, "histograms.root"))
            print(f"{prediction}_OS-SS_lep_{lep}_Dplus_{obs['prior_var']}")
            h = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_{obs['prior_var']}")
            if prediction == "MGFxFx_WplusD":
                gr_qcd = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_{obs['prior_var']}_ratio_mg_fxfx_theory_qcd")
                gr_pdf = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_{obs['prior_var']}_ratio_mg_fxfx_theory_pdf")
            elif prediction == "MGPy8EG_NLO_WplusD":
                gr_qcd = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_{obs['prior_var']}_ratio_mg_nlo_theory_qcd")
                gr_pdf = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_{obs['prior_var']}_ratio_mg_nlo_theory_pdf")
            else:
                gr_qcd = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_{obs['prior_var']}_ratio_sherpa2211_theory_qcd")
                gr_pdf = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_{obs['prior_var']}_ratio_sherpa2211_theory_pdf")
            gr_as = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_{obs['prior_var']}_ratio_sherpa2211_theory_as")
            gr_theory = gr_obs.Clone()
            total_xsec = 0.
            total_xsec_up = 0.
            total_xsec_dn = 0.
            for i in range(h.GetNbinsX()):
                gr_theory.GetY()[i] = h.GetBinContent(i + 1)
                if THEORY_DICT[prediction]["offset"] > 0:
                    xc = gr_theory.GetX()[i] + THEORY_DICT[prediction]["offset"] * gr_theory.GetEXhigh()[i]
                else:
                    xc = gr_theory.GetX()[i] + THEORY_DICT[prediction]["offset"] * gr_theory.GetEXlow()[i]
                xl = obs["bins"][i]
                xh = obs["bins"][i + 1]
                w = xh - xl
                gr_theory.GetX()[i] = xc
                if obs["logx"]:
                    x_err_up = ROOT.TMath.Power(10, ROOT.TMath.Log10(xc) + (ROOT.TMath.Log10(obs["bins"][-1]) - ROOT.TMath.Log10(obs["bins"][0])) / 100.) - xc
                    x_err_dn = xc - ROOT.TMath.Power(10, ROOT.TMath.Log10(xc) - (ROOT.TMath.Log10(obs["bins"][-1]) - ROOT.TMath.Log10(obs["bins"][0])) / 100.)
                else:
                    x_err_up = w / 20.
                    x_err_dn = w / 20.
                gr_theory.GetEXhigh()[i] = x_err_up
                gr_theory.GetEXlow()[i] = x_err_dn
                gr_theory.GetEYhigh()[i] = ((h.GetBinError(i + 1))**2 +
                                            (gr_qcd.GetEYhigh()[i] * gr_theory.GetY()[i])**2 +
                                            (gr_pdf.GetEYhigh()[i] * gr_theory.GetY()[i])**2 +
                                            (gr_as.GetEYhigh()[i] * gr_theory.GetY()[i])**2)**(0.5)
                gr_theory.GetEYlow()[i] = ((h.GetBinError(i + 1))**2 +
                                           (gr_qcd.GetEYlow()[i] * gr_theory.GetY()[i])**2 +
                                           (gr_pdf.GetEYlow()[i] * gr_theory.GetY()[i])**2 +
                                           (gr_as.GetEYlow()[i] * gr_theory.GetY()[i])**2)**(0.5)
                total_xsec += gr_theory.GetY()[i]
                total_xsec_up += gr_theory.GetY()[i] + gr_theory.GetEYhigh()[i]
                total_xsec_dn += gr_theory.GetY()[i] - gr_theory.GetEYlow()[i]

            # style
            gr_theory.SetLineWidth(1)
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
        c1, pad1, pad2, pad3 = createCanvasPads(f"W{lep}_{obs_name}")

        # upper canvas
        pad1.cd()
        if obs["logx"]:
            pad1.SetLogx()
        mg_obs.Draw("a")
        mg_obs.GetXaxis().SetLimits(obs["bins"][0], obs["bins"][-1])
        mg_obs.SetMinimum(1e-3)
        mg_obs.SetMaximum(Y_MAX)

        # ATLAS label
        l1 = ROOT.TLatex()
        l1.SetNDC()
        l1.SetTextFont(73)
        l1.SetTextSize(36)
        l1.DrawLatex(0.19, 0.8 - 0 * 0.08, "ATLAS")
        l2 = ROOT.TLatex()
        l2.SetNDC()
        l2.SetTextFont(43)
        l2.SetTextSize(36)
        l2.DrawLatex(0.305, 0.8 - 0 * 0.08, "Internal")
        l2.DrawLatex(0.19, 0.8 - 1 * 0.08, "#sqrt{s} = 13 TeV, 139 fb^{-1}")
        l2.DrawLatex(0.19, 0.8 - 2 * 0.08, "#it{W}+#it{D}(#rightarrowK#pi#pi), #it{W%s} channel" % ("-" if lep == "minus" else "+"))

        # vertical lines
        lines = []
        for i, x in enumerate(obs["bins"][1:-1]):
            y = 40
            if i > 2:
                y = 20
            line = ROOT.TLine(x, 1e-3, x, y)
            line.SetLineStyle(2)
            line.Draw()
            lines += [line]

        # legend
        leg.Draw()

        ROOT.gPad.RedrawAxis()
        pad2.cd()
        if obs["logx"]:
            pad2.SetLogx()
        # pad2.SetLogy()
        mg_obs_norm.Draw("a")
        mg_obs_norm.GetXaxis().SetLimits(obs["bins"][0], obs["bins"][-1])
        # mg_obs_norm.SetMinimum(0.015)
        mg_obs_norm.SetMinimum(1e-3)
        mg_obs_norm.SetMaximum(0.35)

        # vertical lines
        for x in obs["bins"][1:-1]:
            line = ROOT.TLine(x, 1e-3, x, 0.35)
            line.SetLineStyle(2)
            line.Draw()
            lines += [line]

        ROOT.gPad.RedrawAxis()
        pad3.cd()
        if obs["logx"]:
            pad3.SetLogx()
        pad3.SetGridy()
        mg_obs_ratio.Draw("a")
        mg_obs_ratio.GetXaxis().SetLimits(obs["bins"][0], obs["bins"][-1])
        mg_obs_ratio.SetMinimum(0.83)
        mg_obs_ratio.SetMaximum(1.17)
        line3 = ROOT.TLine(obs["bins"][0], 1, obs["bins"][-1], 1)
        line3.SetLineStyle(2)
        line3.Draw()

        # vertical lines
        for x in obs["bins"][1:-1]:
            line = ROOT.TLine(x, 0.83, x, 1.17)
            line.SetLineStyle(2)
            line.Draw()
            lines += [line]

        ROOT.gPad.RedrawAxis()
        c1.Print(f"fit_results/W{lep}_{obs_name}.pdf")

        trash += [mg_obs]
        trash += [mg_obs_norm]
        trash += [mg_obs_ratio]

    # --------------------------------------------
    # Ladder plots
    # --------------------------------------------
    ROOT.gStyle.SetEndErrorSize(8)
    proxy_pdf = ROOT.TGraph()
    proxy_pdf.SetLineColor(ROOT.kBlack)
    proxy_pdf.SetLineWidth(2)
    proxy_total = ROOT.TGraph()
    proxy_total.SetLineColor(ROOT.kRed)
    proxy_total.SetLineWidth(2)

    # divide by 2 for inclusive W- and W+
    SF = 0.5

    # save for Rc plot
    data = {}

    for lep in ["minus", "plus"]:
        data[lep] = {}

        # data
        if lep == "minus":
            xsec = float(POIs_obs[f"mu_W{lep}_tot"][0]) * priors[f"W{lep}"] * SF
            xsec_err_up = float(POIs_obs[f"mu_W{lep}_tot"][1]) * priors[f"W{lep}"] * SF
            xsec_err_dn = abs(float(POIs_obs[f"mu_W{lep}_tot"][2]) * priors[f"W{lep}"]) * SF
            xsec_err_stat_up = float(POIs_stat[f"mu_W{lep}_tot"][1]) * priors[f"W{lep}"] * SF
            xsec_err_stat_dn = abs(float(POIs_stat[f"mu_W{lep}_tot"][2]) * priors[f"W{lep}"]) * SF
        else:
            xsec = float(POIs_obs[f"expr_mu_W{lep}_tot"][0]) * priors[f"W{lep}"] * SF
            xsec_err_up = float(POIs_obs[f"expr_mu_W{lep}_tot"][1]) * priors[f"W{lep}"] * SF
            xsec_err_dn = abs(float(POIs_obs[f"expr_mu_W{lep}_tot"][1]) * priors[f"W{lep}"]) * SF
            xsec_err_stat_up = float(POIs_stat[f"expr_mu_W{lep}_tot"][1]) * priors[f"W{lep}"] * SF
            xsec_err_stat_dn = abs(float(POIs_stat[f"expr_mu_W{lep}_tot"][1]) * priors[f"W{lep}"]) * SF

        gr = ROOT.TGraph()
        gr.SetPoint(0, xsec, 1.0)
        gr.SetPoint(1, xsec, -1.0)
        gr.SetLineColor(ROOT.kBlack)
        gr.SetLineStyle(2)
        gr.SetLineWidth(2)

        gr_tot = ROOT.TGraphAsymmErrors()
        gr_tot.SetPoint(0, xsec, 0.0)
        gr_tot.SetPointError(0, xsec_err_up, xsec_err_dn, 1.0, 1.0)
        gr_tot.SetFillColor(ROOT.kGray + 1)
        gr_tot.SetLineWidth(0)

        gr_stat = ROOT.TGraphAsymmErrors()
        gr_stat.SetPoint(0, xsec, 0.0)
        gr_stat.SetPointError(0, xsec_err_stat_up, xsec_err_stat_dn, 1.0, 1.0)
        gr_stat.SetFillColor(ROOT.kGray + 2)
        gr_stat.SetLineWidth(0)

        # multi graph
        mg = ROOT.TMultiGraph()
        mg.Add(gr_tot, "e2")
        mg.Add(gr_stat, "e2")
        mg.Add(gr, "l")

        # legend
        N = 5 + len(THEORY_DICT)
        leg = ROOT.TLegend(0.42, 0.80 - N * (0.045), 0.98, 0.80)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(42)
        leg.SetTextFont(43)
        leg.AddEntry(gr, "Data", "l")
        leg.AddEntry(gr_stat, "Stat. Unc.", "f")
        leg.AddEntry(gr_tot, "Syst. #oplus Stat.", "f")

        # theory
        for k, prediction in enumerate(THEORY):
            h = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_{obs['prior_var']}")
            if prediction == "MGFxFx_WplusD":
                gr_qcd = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_{obs['prior_var']}_ratio_mg_fxfx_theory_qcd")
                gr_pdf = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_{obs['prior_var']}_ratio_mg_fxfx_theory_pdf")
            elif prediction == "MGPy8EG_NLO_WplusD":
                gr_qcd = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_{obs['prior_var']}_ratio_mg_nlo_theory_qcd")
                gr_pdf = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_{obs['prior_var']}_ratio_mg_nlo_theory_pdf")
            else:
                gr_qcd = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_{obs['prior_var']}_ratio_sherpa2211_theory_qcd")
                gr_pdf = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_{obs['prior_var']}_ratio_sherpa2211_theory_pdf")
            gr_as = f_theory.Get(f"{prediction}_OS-SS_lep_{lep}_Dplus_{obs['prior_var']}_ratio_sherpa2211_theory_as")
            xsec = 0.
            xsec_up = 0.
            xsec_dn = 0.
            xsec_pdf_up = 0.
            xsec_pdf_dn = 0.
            for i in range(h.GetNbinsX()):
                xsec += h.GetBinContent(i + 1) * SF
                qcd_err_up = gr_qcd.GetEYhigh()[i] * h.GetBinContent(i + 1) * SF
                pdf_err_up = gr_pdf.GetEYhigh()[i] * h.GetBinContent(i + 1) * SF
                as_err_up = gr_as.GetEYhigh()[i] * h.GetBinContent(i + 1) * SF
                qcd_err_dn = gr_qcd.GetEYlow()[i] * h.GetBinContent(i + 1) * SF
                pdf_err_dn = gr_pdf.GetEYlow()[i] * h.GetBinContent(i + 1) * SF
                as_err_dn = gr_as.GetEYlow()[i] * h.GetBinContent(i + 1) * SF
                total_err_up = (qcd_err_up**2 + pdf_err_up**2 + as_err_up**2)**(0.5)
                total_err_dn = (qcd_err_dn**2 + pdf_err_dn**2 + as_err_dn**2)**(0.5)
                xsec_up += total_err_up
                xsec_dn += total_err_dn
                xsec_pdf_up += pdf_err_up
                xsec_pdf_dn += pdf_err_dn

            # save
            data[lep][prediction] = [xsec, xsec_up, xsec_dn, xsec_pdf_up, xsec_pdf_dn]

            # style
            gr_theory_marker = ROOT.TGraph()
            gr_theory_marker.SetPoint(0, xsec, 1 - 0.2 * (k + 1))
            gr_theory_marker.SetMarkerSize(3)
            gr_theory_marker.SetMarkerColor(THEORY_DICT[prediction]["lineColor"])
            gr_theory_marker.SetMarkerStyle(THEORY_DICT[prediction]["markerStyleFull"])

            gr_theory = ROOT.TGraphAsymmErrors()
            gr_theory.SetPoint(0, xsec, 1 - 0.2 * (k + 1))
            gr_theory.SetPointError(0, xsec_up, xsec_dn, 0.0, 0.0)
            gr_theory.SetLineWidth(4)
            gr_theory.SetLineColor(ROOT.kRed)

            gr_theory_pdf = ROOT.TGraphAsymmErrors()
            gr_theory_pdf.SetPoint(0, xsec, 1 - 0.2 * (k + 1))
            gr_theory_pdf.SetPointError(0, xsec_pdf_up, xsec_pdf_dn, 0.0, 0.0)
            gr_theory_pdf.SetLineWidth(4)
            gr_theory_pdf.SetLineColor(ROOT.kBlack)

            mg.Add(gr_theory, "e")
            mg.Add(gr_theory_pdf, "e")
            mg.Add(gr_theory_marker, "p")
            leg.AddEntry(gr_theory_marker, THEORY_DICT[prediction]["legendLabel"], "p")

        leg.AddEntry(proxy_pdf, "PDF Unc.", "l")
        leg.AddEntry(proxy_total, "Total Unc.", "l")

        # canvas
        c = ROOT.TCanvas(f"W{lep}Total", f"W{lep}Total", 1200, 1200)
        c.SetLeftMargin(0.05)
        mg.Draw("a")
        mg.GetYaxis().SetLabelSize(0)
        mg.GetXaxis().SetTitle("#sigma^{tot.}_{fid.} [pb]")
        mg.GetXaxis().SetLimits(30, 100)
        mg.SetMinimum(-1)
        mg.SetMaximum(1)

        # ATLAS label
        l1 = ROOT.TLatex()
        l1.SetTextFont(73)
        l1.SetTextSize(42)
        l1.DrawLatex(60, 0.85, "ATLAS")
        l2 = ROOT.TLatex()
        l2.SetTextFont(43)
        l2.SetTextSize(42)
        l2.DrawLatex(73, 0.85 - 0 * 0.10, "Internal")
        l2.DrawLatex(60, 0.85 - 1 * 0.10, "#sqrt{s} = 13 TeV, 139 fb^{-1}")
        l2.DrawLatex(60, 0.85 - 2 * 0.10, "#it{W}+#it{D}(#rightarrowK#pi#pi), #it{W%s} channel" % ("-" if lep == "minus" else "+"))

        # legend
        leg.Draw()

        # print
        ROOT.gPad.RedrawAxis()
        c.Print(f"fit_results/W{lep}_tot_{obs_name}.pdf")
        trash += [mg]

    # ratio plot
    Rc = float(POIs_obs["mu_Rc"][0])
    Rc_up = float(POIs_obs["mu_Rc"][1])
    Rc_dn = abs(float(POIs_obs["mu_Rc"][2]))
    Rc_stat_up = float(POIs_stat["mu_Rc"][1])
    Rc_stat_dn = abs(float(POIs_stat["mu_Rc"][2]))

    gr = ROOT.TGraph()
    gr.SetPoint(0, Rc, 1.0)
    gr.SetPoint(1, Rc, -1.0)
    gr.SetLineColor(ROOT.kBlack)
    gr.SetLineStyle(2)
    gr.SetLineWidth(2)

    gr_tot = ROOT.TGraphAsymmErrors()
    gr_tot.SetPoint(0, Rc, 0.0)
    gr_tot.SetPointError(0, Rc_up, Rc_dn, 1.0, 1.0)
    gr_tot.SetFillColor(ROOT.kGray + 1)
    gr_tot.SetLineWidth(0)

    gr_stat = ROOT.TGraphAsymmErrors()
    gr_stat.SetPoint(0, Rc, 0.0)
    gr_stat.SetPointError(0, Rc_stat_up, Rc_stat_dn, 1.0, 1.0)
    gr_stat.SetFillColor(ROOT.kGray + 2)
    gr_stat.SetLineWidth(0)

    # multi graph
    mg = ROOT.TMultiGraph()
    mg.Add(gr_tot, "e2")
    mg.Add(gr_stat, "e2")
    mg.Add(gr, "l")

    # theory
    for k, prediction in enumerate(THEORY):
        xsec = data["plus"][prediction][0] / data["minus"][prediction][0]
        xsec_up = (data["plus"][prediction][0] + data["plus"][prediction][1]) / (
            data["minus"][prediction][0] + data["minus"][prediction][1])
        xsec_dn = (data["plus"][prediction][0] - data["plus"][prediction][2]) / (
            data["minus"][prediction][0] - data["minus"][prediction][2])
        xsec_pdf_up = (data["plus"][prediction][0] + data["plus"][prediction][3]) / (
            data["minus"][prediction][0] + data["minus"][prediction][3])
        xsec_pdf_dn = (data["plus"][prediction][0] - data["plus"][prediction][4]) / (
            data["minus"][prediction][0] - data["minus"][prediction][4])

        # style
        gr_theory_marker = ROOT.TGraph()
        gr_theory_marker.SetPoint(0, xsec, 1 - 0.2 * (k + 1))
        gr_theory_marker.SetMarkerSize(3)
        gr_theory_marker.SetMarkerColor(THEORY_DICT[prediction]["lineColor"])
        gr_theory_marker.SetMarkerStyle(THEORY_DICT[prediction]["markerStyleFull"])

        gr_theory = ROOT.TGraphAsymmErrors()
        gr_theory.SetPoint(0, xsec, 1 - 0.2 * (k + 1))
        gr_theory.SetPointError(0, xsec_up - xsec, xsec - xsec_dn, 0.0, 0.0)
        gr_theory.SetLineWidth(4)
        gr_theory.SetLineColor(ROOT.kRed)

        gr_theory_pdf = ROOT.TGraphAsymmErrors()
        gr_theory_pdf.SetPoint(0, xsec, 1 - 0.2 * (k + 1))
        gr_theory_pdf.SetPointError(0, xsec_pdf_up - xsec, xsec - xsec_pdf_dn, 0.0, 0.0)
        gr_theory_pdf.SetLineWidth(4)
        gr_theory_pdf.SetLineColor(ROOT.kBlack)

        mg.Add(gr_theory, "e")
        mg.Add(gr_theory_pdf, "e")
        mg.Add(gr_theory_marker, "p")

    # canvas
    c = ROOT.TCanvas(f"W{lep}Rc", f"W{lep}Rc", 1200, 1200)
    c.SetLeftMargin(0.05)
    mg.Draw("a")
    mg.GetYaxis().SetLabelSize(0)
    mg.GetXaxis().SetTitle("Rc")
    mg.GetXaxis().SetLimits(0.94, 1.2)
    mg.SetMinimum(-1)
    mg.SetMaximum(1)

    # ATLAS label
    l1 = ROOT.TLatex()
    l1.SetTextFont(73)
    l1.SetTextSize(42)
    l1.DrawLatex(0.94 + 30 / 70 * (1.2 - 0.94), 0.85, "ATLAS")
    l2 = ROOT.TLatex()
    l2.SetTextFont(43)
    l2.SetTextSize(42)
    l2.DrawLatex(0.94 + 42 / 70 * (1.2 - 0.94), 0.85 - 0 * 0.10, "Internal")
    l2.DrawLatex(0.94 + 30 / 70 * (1.2 - 0.94), 0.85 - 1 * 0.10, "#sqrt{s} = 13 TeV, 139 fb^{-1}")
    l2.DrawLatex(0.94 + 30 / 70 * (1.2 - 0.94), 0.85 - 2 * 0.10, "#it{W}+#it{D}(#rightarrowK#pi#pi)")

    # legend
    leg.Draw()

    # print
    ROOT.gPad.RedrawAxis()
    c.Print(f"fit_results/Rc_{obs_name}.pdf")

    f.Close()
    f_theory.Close()
