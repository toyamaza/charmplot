#!/usr/bin/env python
from charmplot.common import utils
from charmplot.control import channel
from charmplot.control import sample
from charmplot.control import variable
from ctypes import c_double
import logging
import math
import os
import ROOT
import sys
import yaml

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
sample_MG_FTAG = sample.Sample('FTAG_MG', None, **{'add': ['MadGraph'], 'subtract': [], 'legendLabel': 'Fid. Truth + Reco.', 'lineColor': 'ROOT.kBlue'})
sample_MG_TRUTH = sample.Sample('TRUTH_MG', None, **{'add': ['MadGraph'], 'subtract': [], 'legendLabel': 'Fid. Truth', 'lineColor': 'ROOT.kRed'})
sample_Sherpa = sample.Sample('Sherpa', None, **{'add': ['Sherpa'], 'subtract': [], 'legendLabel': 'Sherpa', 'lineColor': 'ROOT.kBlack'})
sample_MG_Prediction = sample.Sample('MG_Prediction', None, **{'add': ['MadGraph'], 'subtract': [], 'legendLabel': 'MG Truth', 'lineColor': 'ROOT.kBlue'})
sample_MG_FTAG_ALL = sample.Sample('FTAG_MG', None, **{'add': ['MadGraph'], 'subtract': [], 'legendLabel': 'Reco.', 'lineColor': 'ROOT.kRed'})
sample_MG_FTAG_FID = sample.Sample('FTAG_MG', None, **{'add': ['MadGraph'], 'subtract': [], 'legendLabel': 'Fid. Truth + Reco.', 'lineColor': 'ROOT.kBlue'})

# variables to plot
variables = [
    variable.Variable("lep_eta", **{
        "label": "#eta(lep)",
        "x_range": [-3.0, 3.0],
        "ratio_range": [0.89, 1.11],
        "rebin": 5,
    }),
    variable.Variable("lep_pt", **{
        "label": "p_{T}(lep)",
        "unit": "GeV",
        "x_range": [0, 200],
        "ratio_range": [0.89, 1.11],
        "rebin": 5,
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
    variable.Variable("Dmeson_pt", **{
        "label": "p_{T}(D+)",
        "unit": "GeV",
        "x_range": [0, 150],
        "ratio_range": [0.5, 1.5],
        "rebin": 10,
    }),
    variable.Variable("Dmeson_eta", **{
        "label": "#eta(D+)",
        "x_range": [-3.0, 3.0],
        "ratio_range": [0.89, 1.11],
        "rebin": 4,
    }),
]

# rename channel
rename_map = {
    "el_minus" : "W^{-}(e^{-}#bar{#nu})+D^{+}",
    "el_plus" : "W^{+}(e^{+}#nu)+D^{-}",
    "mu_minus" : "W^{-}(#mu^{-}#bar{#nu})+D^{+}",
    "mu_plus" : "W^{+}(#mu^{+}#nu)+D^{-}",
}

def main(options):
    f_d = ROOT.TFile(options.data, "READ")
    f_r = ROOT.TFile(options.mc, "READ")
    f_t = ROOT.TFile(options.truth, "READ")
    logging.info(f"read input files {f_d} {f_r} {f_t}")

    # channels
    channels = options.channels.split(",")

    # TREx output
    trex_yields = {}
    if options.trex:
        for c in channels:
            trex_yields[c] = {'Yield' : 0, 'Error' : 0}
            with open(os.path.join(options.trex, f'WCharm_{c}_all_OS-SS/WCharm_{c}_asimov_OS-SS/Tables/Table_postfit.yaml')) as f:
                table = yaml.load(f, Loader=yaml.FullLoader)
                for x in table:
                    print(x['Region'])
                    for s in x['Samples']:
                        if 'Sample' in s and s['Sample'] == "W+D":
                            print(s['Sample'])
                            print(s['Yield'])
                            print(s['Error'])
                            trex_yields[c]['Yield'] += s['Yield']
                            trex_yields[c]['Error'] = math.sqrt((s['Error'])**2 + (trex_yields[c]['Error'])**2)

    # inclusive histograms
    h_d_inc = ROOT.TH1F("h_d_inc", "h_d_inc", len(channels), 0, len(channels))
    h_r_inc = ROOT.TH1F("h_r_inc", "h_r_inc", len(channels), 0, len(channels))
    h_t_inc = ROOT.TH1F("h_t_inc", "h_t_inc", len(channels), 0, len(channels))
    h_r_all_inc = ROOT.TH1F("h_r_all_inc", "h_r_all_inc", len(channels), 0, len(channels))
    h_r_fid_inc = ROOT.TH1F("h_r_fid_inc", "h_r_fid_inc", len(channels), 0, len(channels))
    for i, c in enumerate(channels):
        h_d_inc.GetXaxis().SetBinLabel(i + 1, rename_map[c])
        if trex_yields:
            h_d_inc.SetBinContent(i + 1, trex_yields[c]['Yield'])
            h_d_inc.SetBinError(i + 1, trex_yields[c]['Error'])
        h_r_inc.GetXaxis().SetBinLabel(i + 1, rename_map[c])
        h_t_inc.GetXaxis().SetBinLabel(i + 1, rename_map[c])
        h_r_all_inc.GetXaxis().SetBinLabel(i + 1, rename_map[c])
        h_r_fid_inc.GetXaxis().SetBinLabel(i + 1, rename_map[c])

    for var in variables:

        for i, c in enumerate(channels):
            reco_string = f"2018_{c}_SR"
            mode_string_OS = f"{options.mode}_OS"
            mode_string_SS = f"{options.mode}_SS"
            logging.info(f"reading histograms")

            # OS names
            name_d_OS = f"{reco_string}_{mode_string_OS}_Matched/{reco_string}_{mode_string_OS}_Matched__{var.name}"
            name_r_OS = f"{reco_string}_fid_{mode_string_OS}_fid_Matched/{reco_string}_fid_{mode_string_OS}_fid_Matched__truth_{var.name}"
            name_t_OS = f"{c}_{mode_string_OS}_fid/{c}_{mode_string_OS}_fid__{var.name.replace('Dmeson', 'D')}"
            name_r_all_OS = f"{reco_string}_{mode_string_OS}_Matched/{reco_string}_{mode_string_OS}_Matched__{var.name}"
            name_r_fid_OS = f"{reco_string}_fid_{mode_string_OS}_fid_Matched/{reco_string}_fid_{mode_string_OS}_fid_Matched__truth_{var.name}"

            # SS names
            name_d_SS = f"{reco_string}_{mode_string_SS}_Matched/{reco_string}_{mode_string_SS}_Matched__{var.name}"
            name_r_SS = f"{reco_string}_fid_{mode_string_SS}_fid_Matched/{reco_string}_fid_{mode_string_SS}_fid_Matched__truth_{var.name}"
            name_t_SS = f"{c}_{mode_string_SS}_fid/{c}_{mode_string_SS}_fid__{var.name.replace('Dmeson', 'D')}"
            name_r_all_SS = f"{reco_string}_{mode_string_SS}_Matched/{reco_string}_{mode_string_SS}_Matched__{var.name}"
            name_r_fid_SS = f"{reco_string}_fid_{mode_string_SS}_fid_Matched/{reco_string}_fid_{mode_string_SS}_fid_Matched__truth_{var.name}"

            # logging
            logging.info(name_d_OS)
            logging.info(name_r_OS)
            logging.info(name_t_OS)
            logging.info(name_r_all_OS)
            logging.info(name_r_fid_OS)

            # OS histograms
            h_d = f_d.Get(name_d_OS).Clone(f"h_d_{reco_string}_OS")
            h_r = f_r.Get(name_r_OS).Clone(f"h_r_{reco_string}_OS")
            h_t = f_t.Get(name_t_OS).Clone(f"h_t_{reco_string}_OS")
            h_r_all = f_r.Get(name_r_all_OS).Clone(f"h_r_all_{reco_string}_OS")
            h_r_fid = f_r.Get(name_r_fid_OS).Clone(f"h_r_all_{reco_string}_OS")

            # SS histograms
            h_d_SS = f_d.Get(name_d_SS).Clone(f"h_d_{reco_string}_SS")
            h_r_SS = f_r.Get(name_r_SS).Clone(f"h_r_{reco_string}_SS")
            h_t_SS = f_t.Get(name_t_SS).Clone(f"h_t_{reco_string}_SS")
            h_r_all_SS = f_r.Get(name_r_all_SS).Clone(f"h_r_all_{reco_string}_SS")
            h_r_fid_SS = f_r.Get(name_r_fid_SS).Clone(f"h_r_all_{reco_string}_SS")

            # Subtract
            h_d.Add(h_d_SS, -1.0)
            h_r.Add(h_r_SS, -1.0)
            h_t.Add(h_t_SS, -1.0)
            h_r_all.Add(h_r_all_SS, -1.0)
            h_r_fid.Add(h_r_fid_SS, -1.0)

            if var == variables[0]:
                # pointers to store errors
                err_d = c_double(0)
                err_r = c_double(0)
                err_t = c_double(0)
                err_r_all = c_double(0)
                err_r_fid = c_double(0)
                int_d = h_d.IntegralAndError(0, h_d.GetNbinsX() + 1, err_d)
                int_r = h_r.IntegralAndError(0, h_r.GetNbinsX() + 1, err_r)
                int_t = h_t.IntegralAndError(0, h_t.GetNbinsX() + 1, err_t)
                int_r_all = h_r_all.IntegralAndError(0, h_r_all.GetNbinsX() + 1, err_r_all)
                int_r_fid = h_r_fid.IntegralAndError(0, h_r_fid.GetNbinsX() + 1, err_r_fid)

                # set inclusive histograms
                if not trex_yields:
                    h_d_inc.SetBinContent(i + 1, int_d)
                h_r_inc.SetBinContent(i + 1, int_r)
                h_t_inc.SetBinContent(i + 1, int_t)
                h_r_all_inc.SetBinContent(i + 1, int_r_all)
                h_r_fid_inc.SetBinContent(i + 1, int_r_fid)
                if not trex_yields:
                    h_d_inc.SetBinError(i + 1, err_d.value)
                h_r_inc.SetBinError(i + 1, err_r.value)
                h_t_inc.SetBinError(i + 1, err_t.value)
                h_r_all_inc.SetBinError(i + 1, err_r_all.value)
                h_r_fid_inc.SetBinError(i + 1, err_r_fid.value)

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

            canv = utils.make_canvas_mc_ratio(mc_map[samples[0]], var, chan, "Efficiency", x=800, y=800, ratio_range=[0, 0.49])

            # configure histograms
            canv.configure_histograms(mc_map, True)

            # top pad
            canv.pad1.cd()
            for s in samples:
                mc_map[s].Draw("hist same")

            # make legend
            canv.make_legend(mc_map, samples, print_yields=False)

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

            canv = utils.make_canvas_mc_ratio(mc_map[samples[0]], var, chan, "FidCorrection", x=800, y=800, ratio_range=[0, 0.99])

            # configure histograms
            canv.configure_histograms(mc_map, True)

            # top pad
            canv.pad1.cd()
            for s in samples:
                mc_map[s].Draw("hist same")

            # make legend
            canv.make_legend(mc_map, samples, print_yields=False)

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
            h_mc_tot.SetLineColor(ROOT.kRed)

            # mc stat error
            gr_mc, gr_mc_stat_err, gr_mc_stat_err_only = utils.make_stat_err_and_nominal(h_mc_tot)
            gr_mc.SetLineColor(ROOT.kRed)
            gr_mc.SetMarkerSize(0.8)
            gr_mc.SetMarkerStyle(24)
            gr_mc.SetMarkerColor(ROOT.kRed)
            gr_mc_stat_err.SetLineColor(ROOT.kRed)
            gr_mc_stat_err.SetFillColor(ROOT.kRed - 10)
            gr_mc_stat_err.SetFillStyle(1001)
            gr_mc_stat_err.SetMarkerSize(0.8)
            gr_mc_stat_err.SetMarkerStyle(24)
            gr_mc_stat_err.SetMarkerColor(ROOT.kRed)
            gr_mc_stat_err_only.SetLineColor(ROOT.kRed)
            gr_mc_stat_err_only.SetFillColor(ROOT.kRed - 10)
            gr_mc_stat_err_only.SetFillStyle(1001)
            gr_mc_stat_err_only.SetMarkerSize(0.8)
            gr_mc_stat_err_only.SetMarkerStyle(24)
            gr_mc_stat_err_only.SetMarkerColor(ROOT.kRed)

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

    # inclusive cross-section
    h_u_inc = h_d_inc.Clone("h_u_inc")
    h_u_inc.Multiply(h_t_inc)
    h_u_inc.Divide(h_r_inc)
    h_u_inc.Multiply(h_r_fid_inc)
    h_u_inc.Divide(h_r_all_inc)
    h_u_inc.Scale(1. / float(options.lumi))

    # prediction
    h_p_inc = h_t_inc.Clone(f"h_p_inc")
    h_p_inc.Scale(1. / float(options.lumi))

    # var and channel object
    var = variable.Variable("Region", **{"label": "Region", "ratio_range": [0.49, 1.51]})
    chan = channel.Channel('inclusive', ['W+jets inclusive', 'inclusive'], "2018", [], [])

    # samples
    samples = [sample_MG_Prediction]

    # mc map
    mc_map = {
        sample_MG_Prediction: h_p_inc,
    }

    # total mc
    hs = utils.make_stack(samples, mc_map)
    h_mc_tot = utils.make_mc_tot(hs, "mc_tot_inclusive")
    h_mc_tot.SetLineColor(ROOT.kRed)

    # mc stat error
    gr_mc, gr_mc_stat_err, gr_mc_stat_err_only = utils.make_stat_err_and_nominal(h_mc_tot)
    gr_mc.SetLineColor(ROOT.kRed)
    gr_mc.SetMarkerSize(0.8)
    gr_mc.SetMarkerStyle(24)
    gr_mc.SetMarkerColor(ROOT.kRed)
    gr_mc_stat_err.SetLineColor(ROOT.kRed)
    gr_mc_stat_err.SetFillColor(ROOT.kRed - 10)
    gr_mc_stat_err.SetFillStyle(1001)
    gr_mc_stat_err.SetMarkerSize(0.8)
    gr_mc_stat_err.SetMarkerStyle(24)
    gr_mc_stat_err.SetMarkerColor(ROOT.kRed)
    gr_mc_stat_err_only.SetLineColor(ROOT.kRed)
    gr_mc_stat_err_only.SetFillColor(ROOT.kRed - 10)
    gr_mc_stat_err_only.SetFillStyle(1001)
    gr_mc_stat_err_only.SetMarkerSize(0.8)
    gr_mc_stat_err_only.SetMarkerStyle(24)
    gr_mc_stat_err_only.SetMarkerColor(ROOT.kRed)

    # data stat error
    gr_data, gr_data_stat_err, gr_data_stat_err_only = utils.make_stat_err_and_nominal(h_u_inc)
    gr_data.SetMarkerSize(0.8)
    gr_data_stat_err.SetFillColor(ROOT.kGray + 3)
    gr_data_stat_err.SetFillStyle(3354)
    gr_data_stat_err.SetMarkerSize(0.8)
    gr_data_stat_err_only.SetFillColor(ROOT.kGray + 3)
    gr_data_stat_err_only.SetFillStyle(3354)
    gr_data_stat_err_only.SetMarkerSize(0.8)
    ROOT.gStyle.SetHatchesSpacing(0.50)

    # canvas
    canv = utils.make_canvas_unfold(h_u_inc, var, chan, x=800, y=800, events="#sigma_{fid}")

    # configure histograms
    canv.configure_histograms(mc_map, h_u_inc)
    h_ratio = utils.make_ratio(h_mc_tot, h_u_inc)
    gr_mc_ratio = gr_mc_stat_err_only.Clone()
    for i in range(0, h_mc_tot.GetNbinsX() + 2):
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
    canv.make_legend(gr_data_stat_err, [gr_mc_stat_err], samples, data_name="Asimov (MG Unf.)")

    # set maximum after creating legend
    canv.set_maximum((h_u_inc, h_mc_tot), var, mc_min=utils.get_mc_min(mc_map, samples))

    # bottom pad
    canv.pad2.cd()
    canv.set_axis_title(canv.proxy_dn, "MC / Data")
    gr_data_stat_err_only.Draw("e2")
    gr_mc_stat_err_only.Draw("e2")
    gr_mc_ratio.Draw("pe")

    # Print out
    canv.print_all(options.output, "inclusive", "integral")


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
                      default="el_plus,el_minus,mu_plus,mu_minus")
    parser.add_option('-l', '--lumi',
                      action="store", dest="lumi",
                      help="luminosity",
                      default=58450.1)
    parser.add_option('--mode',
                      action="store", dest="mode",
                      help="charm meson mode",
                      default="Dplus")
    parser.add_option('--trex',
                      action="store", dest="trex",
                      help="trex output")

    # parse input arguments
    options, args = parser.parse_args()

    # make output folder if not exist
    if not os.path.isdir(options.output):
        os.makedirs(options.output)

    # do the plotting
    main(options)
