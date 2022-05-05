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

# parameters = ROOT.TFile("/global/cfs/cdirs/atlas/wcharm/charmpp_data/systematics/track_eff_systematics.v3.root", "READ")
parameters = ROOT.TFile("/global/cfs/cdirs/atlas/wcharm/charmpp_data/systematics/spg_fit_track_eff_sys.root", "READ")

VARIATIONS = ["nominal", "TRACK_EFF_Overal", "TRACK_EFF_IBL", "TRACK_EFF_PP0", "TRACK_EFF_QGSP", "TOTAL"]


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
    for meson in ["DMinus", "DPlus", "DstarPlus", "DstarMinus"]:

        gr_mean = {var: ROOT.TGraphErrors() for var in VARIATIONS}
        gr_mean_diff = {var: ROOT.TGraphErrors() for var in VARIATIONS}

        meson_lab = {
            "DMinus": "D-",
            "DPlus": "D+",
            "DstarMinus": "D*-",
            "DstarPlus": "D*+",
        }[meson]

        MeV = 1000.0
        inclusive_meson_lab = "Dplus"
        if "star" in meson:
            MeV = 1.0
            inclusive_meson_lab = "Dstar"

        for j, var in enumerate(VARIATIONS):

            if var != "TOTAL":
                gr_mean[var].SetMarkerColor([ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kGray + 1, ROOT.kBlack][j])
                gr_mean[var].SetLineColor([ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kGray + 1, ROOT.kBlack][j])
                gr_mean_diff[var].SetMarkerColor([ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kGray + 1, ROOT.kBlack][j])
                gr_mean_diff[var].SetLineColor([ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kGray + 1, ROOT.kBlack][j])

                gr_mean[var].SetLineWidth(2)
                gr_mean_diff[var].SetLineWidth(2)
            else:
                gr_mean[var].SetLineWidth(0)
                gr_mean[var].SetLineColor(ROOT.kWhite)
                gr_mean[var].SetFillColor(ROOT.kGray + 2)
                gr_mean[var].SetFillStyle(3356)
                gr_mean_diff[var].SetLineWidth(0)
                gr_mean_diff[var].SetLineColor(ROOT.kWhite)
                gr_mean_diff[var].SetFillColor(ROOT.kGray + 2)
                gr_mean_diff[var].SetFillStyle(3356)

            for i in range(1, 6):

                print(f"pars_{meson}_{var}_pt_bin{i}")
                pars_nominal = parameters.Get(f"pars_{meson}_nominal_pt_bin{i}")
                pars = parameters.Get(f"pars_{meson}_{var}_pt_bin{i}")
                pars_total = parameters.Get(f"{inclusive_meson_lab}_TRACK_EFF_TOT_pt_bin{i}_mean_diff")

                if j == 0:
                    gr_mean[var].SetPoint(gr_mean[var].GetN(), i, pars_nominal.GetBinContent(3))
                    gr_mean[var].SetPointError(gr_mean[var].GetN() - 1, 0.5, pars_nominal.GetBinError(3))
                    gr_mean_diff[var].SetPoint(gr_mean_diff[var].GetN(), i, 0)
                    gr_mean_diff[var].SetPointError(gr_mean_diff[var].GetN() - 1, 0.5, 0)
                elif j != len(VARIATIONS) - 1:
                    gr_mean[var].SetPoint(gr_mean[var].GetN(), i + 0.5 - 0.2 * j, pars.GetBinContent(3))
                    gr_mean[var].SetPointError(gr_mean[var].GetN() - 1, 0, pars.GetBinError(3))
                    gr_mean_diff[var].SetPoint(gr_mean_diff[var].GetN(), i + 0.5 - 0.2 * j, MeV * (pars.GetBinContent(3) - pars_nominal.GetBinContent(3)))
                    gr_mean_diff[var].SetPointError(gr_mean_diff[var].GetN() - 1, 0, MeV * (pars.GetBinError(3)**2 + pars_nominal.GetBinError(3)**2)**0.5)
                else:
                    gr_mean[var].SetPoint(gr_mean[var].GetN(), i, pars_nominal.GetBinContent(3))
                    gr_mean[var].SetPointError(gr_mean[var].GetN() - 1, 0.5, pars_total.getVal() / MeV)
                    gr_mean_diff[var].SetPoint(gr_mean_diff[var].GetN(), i, 0)
                    gr_mean_diff[var].SetPointError(gr_mean_diff[var].GetN() - 1, 0.5, pars_total.getVal())

        leg = make_legend(len(VARIATIONS), 0.17, 0.50, 0.70, size=38)
        mg_mean = ROOT.TMultiGraph()
        mg_mean_diff = ROOT.TMultiGraph()
        for var in VARIATIONS:
            if var != "TOTAL":
                leg.AddEntry(gr_mean[var], var, "pe")
            else:
                leg.AddEntry(gr_mean_diff[var], var, "f")

        for var in reversed(["TRACK_EFF_Overal", "TRACK_EFF_IBL", "TRACK_EFF_PP0", "TRACK_EFF_QGSP", "nominal", "TOTAL"]):
            if var != "TOTAL":
                mg_mean.Add(gr_mean[var], "pe")
                mg_mean_diff.Add(gr_mean_diff[var], "pe")
            else:
                mg_mean.Add(gr_mean[var], "e5")
                mg_mean_diff.Add(gr_mean_diff[var], "e5")

        garbage += [mg_mean]
        garbage += [mg_mean_diff]

        c, pad1, pad2 = createCanvasPads(f"mean_{meson}")

        pad1.cd()
        pad1.SetGridy()

        # white box
        if "star" not in meson:
            box = ROOT.TBox(0.4, 1.8705, 3.0, 1.8745)
        else:
            box = ROOT.TBox(0.4, 145.485, 3.0, 145.565)
        box.SetLineWidth(0)
        box.SetFillColor(ROOT.kWhite)
        box.Draw()

        mg_mean.Draw("a")

        atlas_label(["SPG material variations", f"{meson_lab} meson"], "Simulation Internal")
        leg.Draw()
        ROOT.gPad.RedrawAxis()

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

        mg_mean.GetYaxis().SetTitle("mean [MeV]")
        mg_mean.SetMinimum(145.44001)
        mg_mean.SetMaximum(145.57)
        mg_mean_diff.SetMinimum(-0.029)
        mg_mean_diff.SetMaximum(+0.029)
        if "star" not in meson:
            mg_mean.GetYaxis().SetTitle(f"{meson_lab} mean [GeV]")
            mg_mean.SetMinimum(1.868001)
            mg_mean.SetMaximum(1.875)
            mg_mean_diff.SetMinimum(-0.49)
            mg_mean_diff.SetMaximum(+0.49)

        ROOT.gPad.RedrawAxis()
        c.Print(f"spg_variations/{meson}_mean.pdf")

    #
    # resolution
    #
    for meson in ["DMinus", "DPlus", "DstarPlus", "DstarMinus"]:

        gr_sigma = {var: ROOT.TGraphErrors() for var in VARIATIONS}
        gr_sigma_diff = {var: ROOT.TGraphErrors() for var in VARIATIONS}

        meson_lab = {
            "DMinus": "D-",
            "DPlus": "D+",
            "DstarMinus": "D*-",
            "DstarPlus": "D*+",
        }[meson]

        MeV = 1000.0
        inclusive_meson_lab = "Dplus"
        if "star" in meson:
            MeV = 1.0
            inclusive_meson_lab = "Dstar"

        for j, var in enumerate(VARIATIONS):

            if var != "TOTAL":
                gr_sigma[var].SetMarkerColor([ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kGray + 1, ROOT.kBlack][j])
                gr_sigma[var].SetLineColor([ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kGray + 1, ROOT.kBlack][j])
                gr_sigma_diff[var].SetMarkerColor([ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kGray + 1, ROOT.kBlack][j])
                gr_sigma_diff[var].SetLineColor([ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kGray + 1, ROOT.kBlack][j])

                gr_sigma[var].SetLineWidth(2)
                gr_sigma_diff[var].SetLineWidth(2)
            else:
                gr_sigma[var].SetLineWidth(0)
                gr_sigma[var].SetLineColor(ROOT.kWhite)
                gr_sigma[var].SetFillColor(ROOT.kGray + 2)
                gr_sigma[var].SetFillStyle(3356)
                gr_sigma_diff[var].SetLineWidth(0)
                gr_sigma_diff[var].SetLineColor(ROOT.kWhite)
                gr_sigma_diff[var].SetFillColor(ROOT.kGray + 2)
                gr_sigma_diff[var].SetFillStyle(3356)

            for i in range(1, 6):

                print(f"pars_{meson}_{var}_pt_bin{i}")
                pars_nominal = parameters.Get(f"pars_{meson}_nominal_pt_bin{i}")
                pars = parameters.Get(f"pars_{meson}_{var}_pt_bin{i}")
                pars_total = parameters.Get(f"{inclusive_meson_lab}_TRACK_EFF_TOT_pt_bin{i}_sigma_diff")

                if j == 0:
                    gr_sigma[var].SetPoint(i - 1, i, MeV * pars_nominal.GetBinContent(6))
                    gr_sigma[var].SetPointError(i - 1, 0.5, MeV * pars_nominal.GetBinError(6))
                    gr_sigma_diff[var].SetPoint(i - 1, i, 0)
                    gr_sigma_diff[var].SetPointError(i - 1, 0.5, 0)
                elif j != len(VARIATIONS) - 1:
                    gr_sigma[var].SetPoint(i - 1, i + 0.5 - 0.2 * j, MeV * pars.GetBinContent(6))
                    gr_sigma[var].SetPointError(i - 1, 0, MeV * pars.GetBinError(6))

                    # error
                    absdiff = (abs(pars.GetBinContent(6)**2 - pars_nominal.GetBinContent(6)**2))**0.5
                    sign = (pars.GetBinContent(6) - pars_nominal.GetBinContent(6)) / abs(pars.GetBinContent(6) - pars_nominal.GetBinContent(6))
                    err = (pars.GetBinError(6)**2 * pars.GetBinContent(6)**2 / absdiff**2 +
                           pars_nominal.GetBinError(6)**2 * pars_nominal.GetBinContent(6)**2 / absdiff**2)**0.5
                    gr_sigma_diff[var].SetPoint(i - 1, i + 0.5 - 0.2 * j, MeV * sign * absdiff)
                    gr_sigma_diff[var].SetPointError(i - 1, 0, MeV * err)
                else:
                    print(pars_total.getVal(), f"{inclusive_meson_lab}_TRACK_EFF_TOT_pt_bin{i}_sigma_diff")
                    gr_sigma[var].SetPoint(gr_sigma[var].GetN(), i, MeV * pars_nominal.GetBinContent(6))
                    gr_sigma[var].SetPointError(gr_sigma[var].GetN() - 1, 0.5, ((MeV * pars_nominal.GetBinContent(6))**2 + pars_total.getVal()**2)**0.5 - MeV * pars_nominal.GetBinContent(6))
                    gr_sigma_diff[var].SetPoint(gr_sigma_diff[var].GetN(), i, 0)
                    gr_sigma_diff[var].SetPointError(gr_sigma_diff[var].GetN() - 1, 0.5, pars_total.getVal())

        leg = make_legend(len(VARIATIONS), 0.17, 0.50, 0.70, size=38)
        mg_sigma = ROOT.TMultiGraph()
        mg_sigma_diff = ROOT.TMultiGraph()
        for var in VARIATIONS:
            if var != "TOTAL":
                leg.AddEntry(gr_sigma[var], var, "pe")
            else:
                leg.AddEntry(gr_sigma_diff[var], var, "f")

        for var in reversed(["TRACK_EFF_Overal", "TRACK_EFF_IBL", "TRACK_EFF_PP0", "TRACK_EFF_QGSP", "nominal", "TOTAL"]):
            if var != "TOTAL":
                mg_sigma.Add(gr_sigma[var], "pe")
                mg_sigma_diff.Add(gr_sigma_diff[var], "pe")
            else:
                mg_sigma.Add(gr_sigma[var], "e5")
                mg_sigma_diff.Add(gr_sigma_diff[var], "e5")

        garbage += [mg_sigma]
        garbage += [mg_sigma_diff]

        c, pad1, pad2 = createCanvasPads(f"width_{meson}")

        pad1.cd()
        pad1.SetGridy()
        # white box
        if "star" not in meson:
            box = ROOT.TBox(0.4, 21, 3.0, 28.5)
        else:
            box = ROOT.TBox(0.4, 0.92, 3.0, 1.25)
        box.SetLineWidth(0)
        box.SetFillColor(ROOT.kWhite)
        box.Draw()

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
        mg_sigma.SetMinimum(0.7001)
        mg_sigma.SetMaximum(1.3)
        mg_sigma_diff.SetMinimum(-0.49)
        mg_sigma_diff.SetMaximum(+0.49)
        if "star" not in meson:
            mg_sigma.SetMinimum(16.001)
            mg_sigma.SetMaximum(29)
            mg_sigma_diff.SetMinimum(-5.9)
            mg_sigma_diff.SetMaximum(+5.9)

        ROOT.gPad.RedrawAxis()
        c.Print(f"spg_variations/{meson}_width.pdf")


if __name__ == "__main__":

    # make output dir
    if not os.path.isdir("spg_variations"):
        os.makedirs("spg_variations")

    # run
    main()
