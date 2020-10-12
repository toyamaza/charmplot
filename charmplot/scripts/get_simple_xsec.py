#!/usr/bin/env python
from charmplot.common import utils
from charmplot.common import www
from charmplot.control import channel
from charmplot.control import sample
from charmplot.control import variable
from ctypes import c_double
import logging
import math
import os
import ROOT
import sys

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.SetAtlasStyle()

# logging
root = logging.getLogger()
root.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
root.addHandler(handler)


# horizontal error bars for histograms
ROOT.gStyle.SetErrorX(0.5)

# proxy samples
sample_MG_FTAG = sample.Sample('FTAG_MG', None, **{'add': ['MadGraph'], 'subtract': [], 'legendLabel': 'MG FTAG4', 'lineColor': 'ROOT.kBlue'})
sample_MG_TRUTH = sample.Sample('TRUTH_MG', None, **{'add': ['MadGraph'], 'subtract': [], 'legendLabel': 'MG TRUTH1', 'lineColor': 'ROOT.kRed'})
sample_Sherpa = sample.Sample('Sherpa', None, **{'add': ['Sherpa'], 'subtract': [], 'legendLabel': 'Sherpa', 'lineColor': 'ROOT.kBlack'})
sample_MG_Prediction = sample.Sample('MG_Prediction', None, **{'add': ['MadGraph'], 'subtract': [], 'legendLabel': 'MG TRUTH1', 'lineColor': 'ROOT.kBlue'})
sample_MG_FTAG_ALL = sample.Sample('FTAG_MG', None, **{'add': ['MadGraph'], 'subtract': [], 'legendLabel': 'MG All Reco', 'lineColor': 'ROOT.kRed'})
sample_MG_FTAG_FID = sample.Sample('FTAG_MG', None, **{'add': ['MadGraph'], 'subtract': [], 'legendLabel': 'MG Fid Reco', 'lineColor': 'ROOT.kBlue'})

# variables to plot
variables = [
    variable.Variable("lep_eta", **{
        "label": "#eta(lep)",
        "x_range": [-3.0, 3.0],
        "ratio_range": [0.89, 1.11],
        "rebin": 4,
    }),
    variable.Variable("lep_pt", **{
        "label": "p_{T}(lep)",
        "unit": "GeV",
        "x_range": [0, 200],
        "ratio_range": [0.89, 1.11],
        "rebin": 4,
    }),
    variable.Variable("met_mt", **{
        "label": "m_{T}",
        "unit": "GeV",
        "x_range": [0, 200],
        "ratio_range": [0.89, 1.11],
        "rebin": 4,
    }),
    variable.Variable("met_met", **{
        "label": "MET",
        "unit": "GeV",
        "x_range": [0, 200],
        "ratio_range": [0.89, 1.11],
        "rebin": 4,
    }),
]


def main(options):
    f_d = ROOT.TFile(options.data, "READ")
    f_r = ROOT.TFile(options.mc, "READ")
    f_t = ROOT.TFile(options.truth, "READ")
    logging.info(f"read input files {f_d} {f_r} {f_t}")

    # channels
    channels = options.channels.split(",")

    # inclusive histograms
    h_inclusive_d = ROOT.TH1F("h_inclusive_d", "h_inclusive_d", len(channels), 0, len(channels))
    h_inclusive_r = ROOT.TH1F("h_inclusive_r", "h_inclusive_r", len(channels), 0, len(channels))
    h_inclusive_t = ROOT.TH1F("h_inclusive_t", "h_inclusive_t", len(channels), 0, len(channels))
    for i, c in enumerate(channels):
        h_inclusive_d.GetXaxis().SetBinLabel(i + 1, c)
        h_inclusive_r.GetXaxis().SetBinLabel(i + 1, c)
        h_inclusive_t.GetXaxis().SetBinLabel(i + 1, c)

    for var in variables:

        for i, c in enumerate(channels):
            reco_string = f"2018_{c}_SR"
            h_d = f_d.Get(f"{reco_string}/{reco_string}__{var.name}").Clone(f"h_d_{reco_string}")
            h_r = f_r.Get(f"{reco_string}_fid/{reco_string}_fid__truth_{var.name}").Clone(f"h_r_{reco_string}")
            h_t = f_t.Get(f"{c}_fid/{c}_fid__{var.name}").Clone(f"h_t_{reco_string}")
            h_r_all = f_r.Get(f"{reco_string}/{reco_string}__{var.name}").Clone(f"h_r_all_{reco_string}")
            h_r_fid = f_r.Get(f"{reco_string}_fid/{reco_string}_fid__{var.name}").Clone(f"h_r_all_{reco_string}")
            logging.info(f"read histograms {h_d} {h_r} {h_t}")

            if var == variables[0]:
                # pointers to store errors
                err_d = c_double(0)
                err_m = c_double(0)
                err_t = c_double(0)
                int_d = h_d.IntegralAndError(0, h_d.GetNbinsX() + 1, err_d)
                int_m = h_r.IntegralAndError(0, h_r.GetNbinsX() + 1, err_m)
                int_t = h_t.IntegralAndError(0, h_t.GetNbinsX() + 1, err_t)

                # set inclusive histograms
                h_inclusive_d.SetBinContent(i + 1, int_d)
                h_inclusive_r.SetBinContent(i + 1, int_m)
                h_inclusive_t.SetBinContent(i + 1, int_t)
                h_inclusive_d.SetBinError(i + 1, err_d)
                h_inclusive_r.SetBinError(i + 1, err_m)
                h_inclusive_t.SetBinError(i + 1, err_t)

            # rebin
            h_d.Rebin(var.rebin)
            h_r.Rebin(var.rebin)
            h_t.Rebin(var.rebin)
            h_r_all.Rebin(var.rebin)
            h_r_fid.Rebin(var.rebin)

            # channel
            chan = channel.Channel(c, ['W+jets inclusive', c], "2018", [], [])

            # canvas: efficiency
            # -------------------------------

            # samples
            samples = [sample_MG_FTAG, sample_MG_TRUTH]

            # mc map
            mc_map = {
                sample_MG_FTAG: h_r,
                sample_MG_TRUTH: h_t,
            }

            canv = utils.make_canvas_mc_ratio(mc_map[samples[0]], var, chan, "Efficiency", x=800, y=800, ratio_range=[0, 2.0])

            # configure histograms
            canv.configure_histograms(mc_map, True)

            # top pad
            canv.pad1.cd()
            for s in samples:
                mc_map[s].Draw("hist same")

            # make legend
            canv.make_legend(mc_map, samples)

            # set maximum after creating legend
            canv.set_maximum([mc_map[s] for s in samples], var, mc_map[samples[0]])

            # bottom pad
            canv.pad2.cd()

            # ratio histograms
            ratios = []
            for i in range(0, len(samples), 2):
                h = mc_map[samples[i]].Clone(f"{mc_map[samples[i]].GetName()}_ratio")
                h.Divide(mc_map[samples[i + 1]])
                ratios += [h]
                h.Draw("hist same")

            # Print out
            canv.print(f"{options.output}/eff_{c}_{var.name}.pdf")

            # canvas: fiducial correction
            # -------------------------------

            # samples
            samples = [sample_MG_FTAG_FID, sample_MG_FTAG_ALL]

            # mc map
            mc_map = {
                sample_MG_FTAG_FID: h_r_fid,
                sample_MG_FTAG_ALL: h_r_all,
            }

            canv = utils.make_canvas_mc_ratio(mc_map[samples[0]], var, chan, "FidCorrection", x=800, y=800, ratio_range=[0, 2.0])

            # configure histograms
            canv.configure_histograms(mc_map, True)

            # top pad
            canv.pad1.cd()
            for s in samples:
                mc_map[s].Draw("hist same")

            # make legend
            canv.make_legend(mc_map, samples)

            # set maximum after creating legend
            canv.set_maximum([mc_map[s] for s in samples], var, mc_map[samples[0]])

            # bottom pad
            canv.pad2.cd()

            # ratio histograms
            ratios = []
            for i in range(0, len(samples), 2):
                h = mc_map[samples[i]].Clone(f"{mc_map[samples[i]].GetName()}_ratio")
                h.Divide(mc_map[samples[i + 1]])
                ratios += [h]
                h.Draw("hist same")

            # Print out
            canv.print(f"{options.output}/corr_{c}_{var.name}.pdf")

            # canvas: unfolded data
            # -------------------------------
            h_u = h_d.Clone(f"h_u_{reco_string}")
            h_u.Multiply(h_t)
            h_u.Divide(h_r)
            h_u.Multiply(h_r_fid)
            h_u.Divide(h_r_all)
            h_u.Scale(1. / float(options.lumi))

            # prediction
            h_p = h_t.Clone(f"h_p_{reco_string}")
            h_p.Scale(1. / float(options.lumi))

            # samples
            samples = [sample_MG_Prediction]

            # mc map
            mc_map = {
                sample_MG_Prediction: h_p,
            }

            # total mc
            hs = utils.make_stack(samples, mc_map)
            h_mc_tot = utils.make_mc_tot(hs, f"{c}_{var.name}_mc_tot")
            h_mc_tot.SetLineColor(ROOT.kBlue)

            # mc stat error
            gr_mc, gr_mc_stat_err, gr_mc_stat_err_only = utils.make_stat_err_and_nominal(h_mc_tot)
            gr_mc.SetLineColor(ROOT.kBlue)
            gr_mc.SetMarkerSize(0.8)
            gr_mc.SetMarkerStyle(24)
            gr_mc.SetMarkerColor(ROOT.kBlue)
            gr_mc_stat_err.SetLineColor(ROOT.kBlue)
            gr_mc_stat_err.SetFillColor(ROOT.kBlue - 10)
            gr_mc_stat_err.SetFillStyle(1001)
            gr_mc_stat_err.SetMarkerSize(0.8)
            gr_mc_stat_err.SetMarkerStyle(24)
            gr_mc_stat_err.SetMarkerColor(ROOT.kBlue)
            gr_mc_stat_err_only.SetLineColor(ROOT.kBlue)
            gr_mc_stat_err_only.SetFillColor(ROOT.kBlue - 10)
            gr_mc_stat_err_only.SetFillStyle(1001)
            gr_mc_stat_err_only.SetMarkerSize(0.8)
            gr_mc_stat_err_only.SetMarkerStyle(24)
            gr_mc_stat_err_only.SetMarkerColor(ROOT.kBlue)

            # data stat error
            gr_data, gr_data_stat_err, gr_data_stat_err_only = utils.make_stat_err_and_nominal(h_u)
            gr_data.SetMarkerSize(0.8)
            gr_data_stat_err.SetFillColor(ROOT.kGray + 3)
            gr_data_stat_err.SetFillStyle(3354)
            gr_data_stat_err.SetMarkerSize(0.8)
            gr_data_stat_err_only.SetFillColor(ROOT.kGray + 3)
            gr_data_stat_err_only.SetFillStyle(3354)
            gr_data_stat_err_only.SetMarkerSize(0.8)
            ROOT.gStyle.SetHatchesSpacing(0.50)

            # canvas
            canv = utils.make_canvas_unfold(h_u, var, chan, x=800, y=800, events=f"d#sigma / d{var.label}")

            # configure histograms
            canv.configure_histograms(mc_map, h_u)
            h_ratio = utils.make_ratio(h_mc_tot, h_u)
            gr_mc_ratio = gr_mc_stat_err_only.Clone()
            for i in range(0, h_ratio.GetNbinsX() + 2):
                gr_mc_stat_err_only.SetPoint(i, gr_mc_stat_err_only.GetX()[i], h_ratio.GetBinContent(i))
                gr_mc_ratio.SetPoint(i, gr_mc_ratio.GetX()[i], h_ratio.GetBinContent(i))
                gr_mc_ratio.SetPointError(i, h_ratio.GetBinWidth(i) * 0.5, h_ratio.GetBinWidth(i) * 0.5, 0, 0)

            # top pad
            canv.pad1.cd()
            gr_mc_stat_err.Draw("e2")
            gr_mc.Draw("pe")
            gr_data_stat_err.Draw("e2")
            gr_data.Draw("pe")

            # make legend
            canv.make_legend(gr_data_stat_err, [gr_mc_stat_err], samples, data_name="MadGraph (MG Unf.)")

            # set maximum after creating legend
            canv.set_maximum((h_u, h_mc_tot), var, mc_min=utils.get_mc_min(mc_map, samples))

            # bottom pad
            canv.pad2.cd()
            canv.set_axis_title(canv.proxy_dn, "MC / Data")
            gr_data_stat_err_only.Draw("e2")
            gr_mc_stat_err_only.Draw("e2")
            gr_mc_ratio.Draw("pe")

            # Print out
            canv.print_all(options.output, c, var.name)

    # inclusive efficiency
    h_eff = h_inclusive_r.Clone("h_eff")
    h_eff.Divide(h_inclusive_t)

    # inclusive cross-section
    h_xsec = h_inclusive_d.Clone("h_xsec")
    h_xsec.Multiply(h_inclusive_t)
    h_xsec.Divide(h_inclusive_r)
    h_xsec.Scale(1. / float(options.lumi))
    h_xsec.Scale(1. / 1000)  # fb

    # theory cross-section
    gr_xsec = ROOT.TGraphErrors()
    for i in range(h_xsec.GetNbinsX()):
        if i % 2 == 0:
            gr_xsec.SetPoint(i, i + 0.5, 11330 / 1000.)
            gr_xsec.SetPointError(i, 0.5, 300 / 1000.)
        else:
            gr_xsec.SetPoint(i, i + 0.5, 8370 / 1000.)
            gr_xsec.SetPointError(i, 0.5, 220 / 1000.)
    gr_xsec.SetFillColor(ROOT.kGreen + 1)
    gr_xsec.SetLineColor(ROOT.kBlack)

    # plot efficiency
    canv_eff = ROOT.TCanvas("canv_eff", "canv_eff", 800, 800)
    h_eff.Draw("pe text")
    h_eff.SetMinimum(0)
    h_eff.SetMaximum(0.3)
    h_eff.GetYaxis().SetTitle("Total Efficiency")
    h_eff.GetXaxis().SetTitle("Channel")
    ROOT.gPad.RedrawAxis()
    canv_eff.Print(f"{options.output}/inclusive_eff.pdf")

    # plot cross-section
    canv_eff = ROOT.TCanvas("canv_xsec", "canv_xsec", 800, 800)
    h_xsec.Draw("pe")
    gr_xsec.Draw("e2")
    h_xsec.Draw("pe text same")
    h_xsec.SetMinimum(0)
    h_xsec.SetMaximum(20)
    h_xsec.GetYaxis().SetTitle("Total Cross-Section [fb]")
    h_xsec.GetXaxis().SetTitle("Channel")
    ROOT.gPad.RedrawAxis()
    canv_eff.Print(f"{options.output}/inclusive_xsec.pdf")


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-d', '--data',
                      action="store", dest="data",
                      help="data file name")
    parser.add_option('-m', '--mc',
                      action="store", dest="mc",
                      help="mc file name")
    parser.add_option('-t', '--truth',
                      action="store", dest="truth",
                      help="truth file name")
    parser.add_option('-o', '--output',
                      action="store", dest="output",
                      help="output folder name")
    parser.add_option('-c', '--channels',
                      action="store", dest="channels",
                      help="list of channels",
                      default="el_plus,el_minus")
    parser.add_option('-l', '--lumi',
                      action="store", dest="lumi",
                      help="luminosity",
                      default=58450.1)

    # parse input arguments
    options, args = parser.parse_args()

    # make output folder if not exist
    if not os.path.isdir(options.output):
        os.makedirs(options.output)

    # do the plotting
    main(options)
