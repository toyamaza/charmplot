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

FOLDER_RECO = "/pscratch/sd/m/mmuskinj/run/charmpp/v15/validation_2"
FOLDER_TRUTH = "/pscratch/sd/m/mmuskinj/run/charmpp/v15/validation_truth_2"

# COLZ style
ROOT.gStyle.SetPaintTextFormat("1.1e")
ROOT.gStyle.SetPalette(ROOT.kCherry)
ROOT.TColor.InvertPalette()

# LUMI = 58450.1
LUMI = 3244.54 + 33402.2
LUMI_RUN2 = 58450.1 + 44630.6 + 3244.54 + 33402.2

# make output directory
if not os.path.exists("transfer_matrix"):
    os.makedirs("transfer_matrix")

# samples
samples = {
    "Sherpa2211_WplusD_Matched": {
        "color": ROOT.kBlue,
        "reco": "wplusd_sherpa",
        "leg_name": "Sh.2.11 W+D",
    },
    "MGPy8EG_FxFx_WplusD_Matched": {
        "color": ROOT.kRed,
        "reco": "wplusd_madgraph",
        "leg_name": "MG (FxFx) W+D",
    }
}

# variables
vars = {
    "transfer_matrix_track_jet_10_zt": {
        "reco": "Dmeson_transfer_matrix_track_jet_10_zt",
        "truth": "D_jet_10_zt",
        "title": "z_{T}^{#DeltaR=1.0}",
        "xrange": (0.12, 1.2),
        "yrange": (0.12, 1.2),
    },
    "transfer_matrix_track_jet_8_zt": {
        "reco": "Dmeson_transfer_matrix_track_jet_8_zt",
        "truth": "D_jet_8_zt",
        "title": "z_{T}^{#DeltaR=0.8}",
        "xrange": (0.12, 1.2),
        "yrange": (0.12, 1.2),
    },
    "transfer_matrix_track_jet_6_zt": {
        "reco": "Dmeson_transfer_matrix_track_jet_6_zt",
        "truth": "D_jet_6_zt",
        "title": "z_{T}^{#DeltaR=0.6}",
        "xrange": (0.12, 1.2),
        "yrange": (0.12, 1.2),
    },
    "transfer_matrix_track_jet_10_pt": {
        "reco": "Dmeson_transfer_matrix_track_jet_10_pt",
        "truth": "D_jet_10_pt",
        "title": "p_{T}^{#DeltaR=1.0}",
        "xrange": (5, 145),
        "yrange": (5, 145),
    },
    "transfer_matrix_track_jet_8_pt": {
        "reco": "Dmeson_transfer_matrix_track_jet_8_pt",
        "truth": "D_jet_8_pt",
        "title": "p_{T}^{#DeltaR=0.8}",
        "xrange": (5, 145),
        "yrange": (5, 145),
    },
    "transfer_matrix_track_jet_6_pt": {
        "reco": "Dmeson_transfer_matrix_track_jet_6_pt",
        "truth": "D_jet_6_pt",
        "title": "p_{T}^{#DeltaR=0.6}",
        "xrange": (5, 145),
        "yrange": (5, 145),
    },
}


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
    GLOBAL_SF = 1.1
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
    h_dn.GetXaxis().SetLabelSize(h_dn.GetXaxis().GetLabelSize() * SF * GLOBAL_SF)

    # tick marks
    h_up.GetYaxis().SetNdivisions(506)
    h_dn.GetYaxis().SetNdivisions(306)


def main():

    f_truth = ROOT.TFile(os.path.join(FOLDER_TRUTH, "truth/wplusd_truth_analysis/histograms.root"))

    for var_name, var in vars.items():

        eff_map_x = {}
        eff_map_y = {}

        for sample_name, sample in samples.items():

            f_reco = ROOT.TFile(os.path.join(FOLDER_RECO, f"charm_frag/{sample['reco']}/histograms.root"))

            # get the histograms
            h_reco = f_reco.Get(f"{sample_name}_OS-SS_0tag_Dplus_{var['reco']}")
            h_reco.Scale(1. / (LUMI * 0.094))
            h_reco_y = h_reco.ProjectionY()
            h_reco_x = h_reco.ProjectionX()
            h_truth = f_truth.Get(f"{sample_name}_OS-SS_Dplus_Kpipi_{var['truth']}")
            h_truth.Scale(LUMI_RUN2 / LUMI)

            # check for consistency
            assert h_reco.GetNbinsX() == h_truth.GetNbinsX()
            assert h_reco.GetNbinsY() == h_truth.GetNbinsX()

            # rebin
            h_reco.GetXaxis().SetRangeUser(var["xrange"][0], var["xrange"][1])
            h_reco.GetYaxis().SetRangeUser(var["yrange"][0], var["yrange"][1])
            h_reco_y.GetXaxis().SetRangeUser(var["xrange"][0], var["xrange"][1])
            h_reco_x.GetXaxis().SetRangeUser(var["xrange"][0], var["xrange"][1])
            h_truth.GetXaxis().SetRangeUser(var["xrange"][0], var["xrange"][1])
            h_matrix = h_reco.Clone(f"{h_reco.GetName()}_matrix")

            # divide reco by truth to get efficiency
            for i in range(1, h_matrix.GetNbinsX() + 1):
                for j in range(1, h_matrix.GetNbinsY() + 1):
                    if h_truth.GetBinContent(i) == 0 or h_matrix.GetBinContent(i, j) / h_truth.GetBinContent(i) > 1 or h_matrix.GetBinContent(i, j) / h_truth.GetBinContent(i) < 0:
                        h_matrix.SetBinContent(i, j, 0)
                        h_matrix.SetBinError(i, j, 0)
                    else:
                        h_matrix.SetBinContent(i, j, h_matrix.GetBinContent(i, j) / h_truth.GetBinContent(i))

            # # normalized matrix
            # integral = h_matrix.GetSumOfWeights()
            # for i in range(1, h_matrix.GetNbinsX() + 1):
            #     for j in range(1, h_matrix.GetNbinsX() + 1):
            #         z = h_matrix.GetBinContent(i, j)
            #         h_matrix.SetBinContent(i, j, 100 * z / integral)

            # Make the 2D plot
            canv = ROOT.TCanvas("_".join([var_name, sample_name]), "_".join([var_name, sample_name]), 1000, 800)
            canv.SetRightMargin(0.18)
            canv.SetLeftMargin(0.18)
            canv.SetBottomMargin(0.20)
            h_matrix.GetZaxis().SetTitle("Normalized Response Matrix [%]")
            h_matrix.GetXaxis().SetTitle(f"{var['title']} (reco)")
            h_matrix.GetYaxis().SetTitle(f"{var['title']} (truth)")
            h_matrix.GetZaxis().SetTitleSize(h_matrix.GetZaxis().GetTitleSize() * 0.9)
            h_matrix.GetZaxis().SetLabelSize(h_matrix.GetZaxis().GetLabelSize() * 0.97)
            h_matrix.Scale(100.)
            h_matrix.SetMinimum(0.01)
            h_matrix.Draw("colz")

            # ATLAS label
            ROOT.ATLASLabel(0.20, 0.90, "Simulation Internal", 1, 0.04)
            ROOT.myText(0.20, 0.85, 1, "#sqrt{s} = 13 TeV", 0.04)
            ROOT.myText(0.20, 0.80, 1, "W+D(#rightarrowK#pi#pi)", 0.04)
            ROOT.myText(0.20, 0.75, 1, sample_name, 0.04)

            # save
            ROOT.gPad.RedrawAxis()
            canv.Print(f"transfer_matrix/{var_name}_{sample_name}.pdf")

            # Make the 1D plots
            for i in range(0, h_reco_y.GetNbinsX() + 2):
                if (h_reco_y.GetBinContent(i) < 0 or h_truth.GetBinContent(i) < 0) or h_truth.GetBinContent(i) == 0 or (h_reco_y.GetBinContent(i) / h_truth.GetBinContent(i) > 1):
                    h_reco_y.SetBinContent(i, 1e-6)
                    h_reco_y.SetBinError(i, 0)
                    h_truth.SetBinContent(i, 1e6)
                    h_truth.SetBinError(i, 0)

            for i in range(0, h_reco_x.GetNbinsX() + 2):
                if (h_reco_x.GetBinContent(i) < 0 or h_truth.GetBinContent(i) < 0) or h_truth.GetBinContent(i) == 0 or (h_reco_x.GetBinContent(i) / h_truth.GetBinContent(i) > 1):
                    h_reco_x.SetBinContent(i, 1e-6)
                    h_reco_x.SetBinError(i, 0)
                    h_truth.SetBinContent(i, 1e6)
                    h_truth.SetBinError(i, 0)

            # projection y
            eff_y = ROOT.TEfficiency(h_reco_y, h_truth)
            eff_y_gr = eff_y.CreateGraph()
            assert h_reco_y.GetNbinsX() == h_truth.GetNbinsX() and h_reco_y.GetNbinsX() == eff_y_gr.GetN(), f"{h_reco_y.GetNbinsX()} {h_truth.GetNbinsX()} {eff_y_gr.GetN()}"
            eff_y_gr.SetMarkerColor(sample["color"])
            eff_y_gr.SetLineColor(sample["color"])
            eff_map_y[sample_name] = eff_y_gr

            # projection x
            eff_x = ROOT.TEfficiency(h_reco_x, h_truth)
            eff_x_gr = eff_x.CreateGraph()
            assert h_reco_x.GetNbinsX() == h_truth.GetNbinsX() and h_reco_x.GetNbinsX() == eff_x_gr.GetN(), f"{h_reco_x.GetNbinsX()} {h_truth.GetNbinsX()} {eff_x_gr.GetN()}"
            eff_x_gr.SetMarkerColor(sample["color"])
            eff_x_gr.SetLineColor(sample["color"])
            eff_map_x[sample_name] = eff_x_gr

        # MGraphs
        mg_y = ROOT.TMultiGraph()
        mg_y_ratio = ROOT.TMultiGraph()
        mg_x = ROOT.TMultiGraph()
        mg_x_ratio = ROOT.TMultiGraph()

        for sample, gr in eff_map_y.items():
            mg_y.Add(gr, "pe")
            gr_ratio = gr.Clone(f"{gr.GetName()}_{sample}_ratio")
            for i in range(0, gr_ratio.GetN()):
                y_nominal = eff_map_y["Sherpa2211_WplusD_Matched"].GetY()[i]
                if y_nominal == 0:
                    gr_ratio.GetY()[i] = 0.
                    gr_ratio.GetEYhigh()[i] = 0.
                    gr_ratio.GetEYlow()[i] = 0.
                else:
                    gr_ratio.GetY()[i] /= y_nominal
                    gr_ratio.GetEYhigh()[i] /= y_nominal
                    gr_ratio.GetEYlow()[i] /= y_nominal
            mg_y_ratio.Add(gr_ratio, "pe")

        for sample, gr in eff_map_x.items():
            mg_x.Add(gr, "pe")
            gr_ratio = gr.Clone(f"{gr.GetName()}_{sample}_ratio")
            for i in range(0, gr_ratio.GetN()):
                y_nominal = eff_map_x["Sherpa2211_WplusD_Matched"].GetY()[i]
                if y_nominal == 0:
                    gr_ratio.GetY()[i] = 0.
                    gr_ratio.GetEYhigh()[i] = 0.
                    gr_ratio.GetEYlow()[i] = 0.
                else:
                    gr_ratio.GetY()[i] /= y_nominal
                    gr_ratio.GetEYhigh()[i] /= y_nominal
                    gr_ratio.GetEYlow()[i] /= y_nominal
            mg_x_ratio.Add(gr_ratio, "pe")

        # legend
        N = len(samples)
        leg = ROOT.TLegend(0.62, 0.848 - N / 2. * (0.1), 0.92, 0.848)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(42)
        leg.SetTextFont(43)
        for sample_name, sample in samples.items():
            leg.AddEntry(eff_map_y[sample_name], sample["leg_name"], "pe")

        for mg, mg_ratio, x, label in zip([mg_y, mg_x], [mg_y_ratio, mg_x_ratio], ["y", "x"], ["truth", "reco"]):
            # canvas
            canv, pad1, pad2 = createCanvasPads(f"{var_name}_{x}")
            configure_axis(mg, mg_ratio)
            pad1.cd()
            mg.Draw("a")
            mg.GetYaxis().SetTitle("Efficiency")
            mg.SetMinimum(0)
            mg.SetMaximum(max([gr.GetHistogram().GetMaximum() for gr in mg.GetListOfGraphs()]) * 1.2)

            # legend
            leg.Draw()

            # ATLAS label
            l1 = ROOT.TLatex()
            l1.SetTextFont(73)
            l1.SetTextSize(42)
            l1.DrawLatexNDC(0.18, 0.82, "ATLAS")
            l2 = ROOT.TLatex()
            l2.SetTextFont(43)
            l2.SetTextSize(42)
            l2.DrawLatexNDC(0.31, 0.82, "Simulation Internal")
            l2.DrawLatexNDC(0.18, 0.82 - 1 * 0.06, "#sqrt{s} = 13 TeV")
            l2.DrawLatexNDC(0.18, 0.82 - 2 * 0.06, "W+D(#rightarrowK#pi#pi)")
            ROOT.gPad.RedrawAxis()

            # pad2
            pad2.cd()
            mg_ratio.Draw("a")
            mg_ratio.GetXaxis().SetTitle(f"{var['title']} ({label})")
            mg_ratio.GetYaxis().SetTitle("Ratio to Sherpa")
            mg_ratio.SetMinimum(0.51)
            mg_ratio.SetMaximum(1.79)
            ROOT.gPad.RedrawAxis()

            # limits
            mg.GetXaxis().SetLimits(var["xrange"][0], var["xrange"][1])
            mg_ratio.GetXaxis().SetLimits(var["xrange"][0], var["xrange"][1])

            # save
            canv.Print(f"transfer_matrix/{var_name}_{x}.pdf")


if __name__ == "__main__":

    # run
    main()
