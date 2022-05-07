#!/usr/bin/env python
import os
import ROOT
dirname = os.path.join(os.path.dirname(__file__), "../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasLabels.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasUtils.C"))
ROOT.SetAtlasStyle()
ROOT.gStyle.SetEndErrorSize(8)

rand = ROOT.TRandom3()

parameters = ROOT.TFile("fits/Dplus_mass_fit.root", "READ")

# VARIATIONS = ["Sherpa2211_WplusD", "MGPy8EG_NLO_WplusD", "MGPy8EG_FxFx_WplusD", "Sherpa2211_Wjets", "MG_Wjets"]
VARIATIONS = ["Sherpa2211_WplusD", "MGPy8EG_NLO_WplusD", "MGPy8EG_FxFx_WplusD", "MG_Wjets"]


def createCanvasPads(name):
    c = ROOT.TCanvas(name, name, 1200, 1200)
    # Upper histogram plot is pad1
    pad1 = ROOT.TPad(f"pad1_{name}", f"pad1_{name}", 0, 0.30, 1, 1.0)
    pad1.SetTopMargin(0.1)
    pad1.SetBottomMargin(0.05 / 0.70)
    pad1.SetFillStyle(4000)
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()
    pad2 = ROOT.TPad(f"pad2_{name}", f"pad2_{name}", 0, 0.00, 1, 0.40)
    pad2.SetTopMargin(0.05 / 0.40)
    pad2.SetBottomMargin(0.10 / 0.40)
    pad2.SetFillStyle(4000)
    pad2.Draw()

    return c, pad1, pad2


def configure_axis(h_up, h_dn):
    GLOBAL_SF = 1.0
    h_up.GetYaxis().SetTitleSize(h_up.GetYaxis().GetTitleSize() * GLOBAL_SF)
    h_up.GetYaxis().SetTitleOffset(h_up.GetYaxis().GetTitleOffset() * (1 / (GLOBAL_SF * 1.1)))
    h_up.GetYaxis().SetLabelSize(h_up.GetYaxis().GetLabelSize() * GLOBAL_SF)
    h_up.GetXaxis().SetLabelSize(0)

    SF = 0.70 / 0.40
    h_dn.GetYaxis().SetTitleSize(h_dn.GetYaxis().GetTitleSize() * SF * GLOBAL_SF)
    h_dn.GetYaxis().SetTitleOffset(h_dn.GetYaxis().GetTitleOffset() * (1 / (GLOBAL_SF * 1.06 * SF)))
    h_dn.GetXaxis().SetTitleSize(h_dn.GetXaxis().GetTitleSize() * SF * GLOBAL_SF)
    h_dn.GetXaxis().SetTitleOffset(h_dn.GetXaxis().GetTitleOffset() * (1 / (GLOBAL_SF * 0.6 * SF)))
    h_dn.GetXaxis().SetLabelOffset(h_dn.GetXaxis().GetLabelOffset() * 4.0)
    h_dn.GetYaxis().SetLabelSize(h_dn.GetYaxis().GetLabelSize() * SF * GLOBAL_SF)
    h_dn.GetXaxis().SetLabelSize(h_dn.GetXaxis().GetLabelSize() * SF * GLOBAL_SF * 1.2)

    # tick marks
    h_up.GetYaxis().SetNdivisions(506)
    h_dn.GetYaxis().SetNdivisions(306)


def atlas_label(labels=[], internal="Internal"):
    l1 = ROOT.TLatex()
    l1.SetTextFont(73)
    l1.SetTextSize(42)
    l1.DrawLatexNDC(0.18, 0.82, "ATLAS")
    l2 = ROOT.TLatex()
    l2.SetTextFont(43)
    l2.SetTextSize(42)
    l2.DrawLatexNDC(0.33, 0.82, "Internal")
    for i, lab in enumerate(labels):
        l2.DrawLatexNDC(0.18, 0.82 - (i + 1) * 0.05, lab)


def make_legend(N, x1=0.58, x2=0.92, y2=0.848, w=0.06, size=42):
    leg = ROOT.TLegend(x1, y2 - N * w, x2, y2)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(size)
    leg.SetTextFont(43)
    return leg


def main():

    garbage = []

    #
    # mean
    #
    for lepton in ["plus", "minus"]:

        gr_mean = {var: ROOT.TGraphErrors() for var in VARIATIONS}
        gr_mean_diff = {var: ROOT.TGraphErrors() for var in VARIATIONS}

        meson_lab = {
            "plus": "D-",
            "minus": "D+",
        }[lepton]

        for j, var in enumerate(VARIATIONS):

            gr_mean[var].SetMarkerColor([ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGray + 2, ROOT.kGray + 2][j])
            gr_mean[var].SetLineColor([ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGray + 2, ROOT.kGray + 2][j])
            gr_mean_diff[var].SetMarkerColor([ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGray + 2, ROOT.kGray + 2][j])
            gr_mean_diff[var].SetLineColor([ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGray + 2, ROOT.kGray + 2][j])
            gr_mean[var].SetLineWidth(2)
            gr_mean_diff[var].SetLineWidth(2)

            for i in range(1, 6):

                print(f"pars_{var}_Matched_truth_pt_bin{i}_OS-SS_lep_{lepton}_Dplus_truth_pt_bin{i}_pt_bin{i}_Dmeson_m_norebin")
                pars_nominal = parameters.Get(
                    f"pars_Sherpa2211_WplusD_Matched_truth_pt_bin{i}_OS-SS_lep_{lepton}_Dplus_truth_pt_bin{i}_pt_bin{i}_Dmeson_m_norebin")
                pars = parameters.Get(f"pars_{var}_Matched_truth_pt_bin{i}_OS-SS_lep_{lepton}_Dplus_truth_pt_bin{i}_pt_bin{i}_Dmeson_m_norebin")

                if j == 0:
                    gr_mean[var].SetPoint(gr_mean[var].GetN(), i, pars_nominal.GetBinContent(3))
                    gr_mean[var].SetPointError(gr_mean[var].GetN() - 1, 0.5, pars_nominal.GetBinError(3))
                    gr_mean_diff[var].SetPoint(gr_mean_diff[var].GetN(), i, 0)
                    gr_mean_diff[var].SetPointError(gr_mean_diff[var].GetN() - 1, 0.5, 0)
                else:
                    gr_mean[var].SetPoint(gr_mean[var].GetN(), i + 0.5 - 0.2 * j, pars.GetBinContent(3))
                    gr_mean[var].SetPointError(gr_mean[var].GetN() - 1, 0.08, pars.GetBinError(3))
                    gr_mean_diff[var].SetPoint(gr_mean_diff[var].GetN(), i + 0.5 - 0.2 * j, 1000 * (pars.GetBinContent(3) - pars_nominal.GetBinContent(3)))
                    gr_mean_diff[var].SetPointError(gr_mean_diff[var].GetN() - 1, 0.08, 1000 * (pars.GetBinError(3)**2 + pars_nominal.GetBinError(3)**2)**0.5)

        leg = make_legend(len(VARIATIONS), 0.17, 0.50, 0.70, size=38)
        mg_mean = ROOT.TMultiGraph()
        mg_mean_diff = ROOT.TMultiGraph()
        for var in VARIATIONS:
            leg.AddEntry(gr_mean[var], var, "pe")

        for var in VARIATIONS:
            mg_mean.Add(gr_mean[var], "pe")
            mg_mean_diff.Add(gr_mean_diff[var], "pe")

        garbage += [mg_mean]
        garbage += [mg_mean_diff]

        c, pad1, pad2 = createCanvasPads(f"mean_{lepton}")

        pad1.cd()
        pad1.SetGridy()

        # white box
        box = ROOT.TBox(0.4, 1.8705, 3.0, 1.8745)
        box.SetLineWidth(0)
        box.SetFillColor(ROOT.kWhite)
        box.Draw()

        mg_mean.Draw("a")

        atlas_label(["W+D(#rightarrowK#pi#pi), OS-SS", f"{meson_lab} meson"], "Simulation Internal")
        leg.Draw()

        pad2.cd()
        pad2.SetGridy()
        mg_mean_diff.Draw("a")

        configure_axis(mg_mean, mg_mean_diff)
        mg_mean.GetXaxis().SetNdivisions(6)
        mg_mean.GetYaxis().SetNdivisions(508)
        mg_mean.GetYaxis().SetTitleOffset(1.2 * mg_mean.GetYaxis().GetTitleOffset())
        mg_mean_diff.GetXaxis().SetTitle("p_{T} bin number")
        mg_mean_diff.GetYaxis().SetTitle("#Delta mean [MeV]")
        mg_mean_diff.GetXaxis().SetNdivisions(6)

        mg_mean.GetYaxis().SetTitle(f"{meson_lab} mean [GeV]")
        mg_mean.SetMinimum(1.86501)
        mg_mean.SetMaximum(1.88)
        mg_mean_diff.SetMinimum(-1.99)
        mg_mean_diff.SetMaximum(+1.99)
        mg_mean.SetMinimum(1.868001)
        mg_mean.SetMaximum(1.875)
        mg_mean_diff.SetMinimum(-1.49)
        mg_mean_diff.SetMaximum(+1.49)

        c.Print(f"wpluds_peak_shape/{lepton}_mean.pdf")

    #
    # resolution
    #
    for lepton in ["plus", "minus"]:

        gr_sigma = {var: ROOT.TGraphErrors() for var in VARIATIONS}
        gr_sigma_diff = {var: ROOT.TGraphErrors() for var in VARIATIONS}

        meson_lab = {
            "plus": "D-",
            "minus": "D+",
        }[lepton]

        for j, var in enumerate(VARIATIONS):

            gr_sigma[var].SetMarkerColor([ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGray + 2, ROOT.kGray + 2][j])
            gr_sigma[var].SetLineColor([ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGray + 2, ROOT.kGray + 2][j])
            gr_sigma_diff[var].SetMarkerColor([ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGray + 2, ROOT.kGray + 2][j])
            gr_sigma_diff[var].SetLineColor([ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGray + 2, ROOT.kGray + 2][j])

            gr_sigma[var].SetLineWidth(2)
            gr_sigma_diff[var].SetLineWidth(2)

            for i in range(1, 6):

                print(f"pars_{var}_Matched_truth_pt_bin{i}_OS-SS_lep_{lepton}_Dplus_truth_pt_bin{i}_pt_bin{i}_Dmeson_m_norebin")
                pars_nominal = parameters.Get(
                    f"pars_Sherpa2211_WplusD_Matched_truth_pt_bin{i}_OS-SS_lep_{lepton}_Dplus_truth_pt_bin{i}_pt_bin{i}_Dmeson_m_norebin")
                pars = parameters.Get(f"pars_{var}_Matched_truth_pt_bin{i}_OS-SS_lep_{lepton}_Dplus_truth_pt_bin{i}_pt_bin{i}_Dmeson_m_norebin")

                if j == 0:
                    gr_sigma[var].SetPoint(i - 1, i, 1000 * pars_nominal.GetBinContent(6))
                    gr_sigma[var].SetPointError(i - 1, 0.5, 1000 * pars_nominal.GetBinError(6))
                    gr_sigma_diff[var].SetPoint(i - 1, i, 0)
                    gr_sigma_diff[var].SetPointError(i - 1, 0.5, 0)
                else:
                    gr_sigma[var].SetPoint(i - 1, i + 0.5 - 0.2 * j, 1000 * pars.GetBinContent(6))
                    gr_sigma[var].SetPointError(i - 1, 0.08, 1000 * pars.GetBinError(6))

                    # error
                    absdiff = (abs(pars.GetBinContent(6)**2 - pars_nominal.GetBinContent(6)**2))**0.5
                    sign = (pars.GetBinContent(6) - pars_nominal.GetBinContent(6)) / abs(pars.GetBinContent(6) - pars_nominal.GetBinContent(6))
                    err = (pars.GetBinError(6)**2 * pars.GetBinContent(6)**2 / absdiff**2 +
                           pars_nominal.GetBinError(6)**2 * pars_nominal.GetBinContent(6)**2 / absdiff**2)**0.5
                    gr_sigma_diff[var].SetPoint(i - 1, i + 0.5 - 0.2 * j, 1000 * sign * absdiff)
                    gr_sigma_diff[var].SetPointError(i - 1, 0.08, 1000 * err)

        leg = make_legend(len(VARIATIONS), 0.17, 0.50, 0.70, size=38)
        mg_sigma = ROOT.TMultiGraph()
        mg_sigma_diff = ROOT.TMultiGraph()
        for var in VARIATIONS:
            leg.AddEntry(gr_sigma[var], var, "pe")

        for var in VARIATIONS:
            mg_sigma.Add(gr_sigma[var], "pe")
            mg_sigma_diff.Add(gr_sigma_diff[var], "pe")

        garbage += [mg_sigma]
        garbage += [mg_sigma_diff]

        c, pad1, pad2 = createCanvasPads(f"width_{lepton}")

        pad1.cd()
        pad1.SetGridy()
        # # white box
        # if "star" not in meson:
        #     box = ROOT.TBox(0.4, 21, 3.0, 28.5)
        # else:
        #     box = ROOT.TBox(0.4, 0.92, 3.0, 1.25)
        # box.SetLineWidth(0)
        # box.SetFillColor(ROOT.kWhite)
        # box.Draw()

        mg_sigma.Draw("a")

        atlas_label(["SPG material variations", f"{meson_lab} meson"], "Simulation Internal")
        leg.Draw()

        pad2.cd()
        pad2.SetGridy()
        mg_sigma_diff.Draw("a")

        configure_axis(mg_sigma, mg_sigma_diff)
        mg_sigma.GetXaxis().SetNdivisions(6)
        mg_sigma.GetYaxis().SetNdivisions(508)
        mg_sigma_diff.GetXaxis().SetTitle("p_{T} bin number")
        mg_sigma_diff.GetYaxis().SetTitle("sqrt(|#sigma_{0}^{2} - #sigma^{2}|) [MeV]")
        mg_sigma_diff.GetXaxis().SetNdivisions(6)

        mg_sigma.GetYaxis().SetTitle(f"{meson_lab} width (#sigma) [MeV]")
        mg_sigma.SetMinimum(16.001)
        mg_sigma.SetMaximum(35)
        mg_sigma_diff.SetMinimum(-8.9)
        mg_sigma_diff.SetMaximum(+8.9)

        c.Print(f"wpluds_peak_shape/{lepton}_width.pdf")


if __name__ == "__main__":

    # make output dir
    if not os.path.isdir("wpluds_peak_shape"):
        os.makedirs("wpluds_peak_shape")

    # run
    main()
