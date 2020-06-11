#!/usr/bin/env python
from charmplot.common import utils
from charmplot.common import www
from charmplot.control import channel
from charmplot.control import sample
from charmplot.control import variable
import logging
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

# sample object
sample_QCDTemplateFit = sample.Sample('Multijet', None, **{'add': ['Multijet'], 'subtract': [], 'legendLabel': 'QCD Fit'})
sample_MatrixMethod = sample.Sample('Multijet', None, **{'add': ['Multijet'], 'subtract': [], 'legendLabel': 'Matrix M.'})
samples = [
    sample_QCDTemplateFit,
    sample_MatrixMethod,
]

# horizontal error bars for histograms
ROOT.gStyle.SetErrorX(0.5)

# variables to plot
variables = [
    variable.Variable("lep_pt", **{
        "label": "p_{T}",
        "unit": "GeV",
        "x_range": [0, 200],
    }),
    variable.Variable("Dmeson_m", **{
        "label": "D^{#pm} m",
        "unit": "GeV",
    }),
    variable.Variable("Dmeson_pt", **{
        "label": 'D^{#pm} p_{T}',
        "unit": "GeV",
    }),
    variable.Variable("met_met", **{
        "label": "MET",
        "unit": "GeV",
        "x_range": [0, 200],
    }),
    variable.Variable("met_mt", **{
        "label": "m_{T}",
        "unit": "GeV",
        "x_range": [0, 200],
    }),
]

# channel object
channels = [
    channel.Channel("OS_2018_el_SR_Dplus", ['W#rightarrowe#nu+D, D#rightarrowK#pi#pi, OS', 'Signal Region'], "2018", [], []),
    channel.Channel("SS_2018_el_SR_Dplus", ['W#rightarrowe#nu+D, D#rightarrowK#pi#pi, SS', 'Signal Region'], "2018", [], []),
    channel.Channel("OS-SS_2018_el_SR_Dplus", ['W#rightarrowe#nu+D, D#rightarrowK#pi#pi, OS-SS', 'Signal Region'], "2018", [], []),
]

# horizontal error bars for histograms
ROOT.gStyle.SetErrorX(0.5)

# input files
f_QCDTemplateFit = ROOT.TFile("/global/u2/m/mmuskinj/work/run/charmpp/v5/wplusd_qcd_pt_2_rw/likelihood_fit/wplusd/histograms.root")
f_MatrixMethod = ROOT.TFile("/global/u2/m/mmuskinj/work/run/charmpp/v5/FakeRate/wplusd/fake_rate/wplusd_madgraph/histograms.root")


def main(options):
    for chan in channels:

        # keep track of first/last plot of each channel
        first_plot = True

        # make channel folder if not exist
        if not os.path.isdir(os.path.join(options.out_name, chan.name)):
            os.makedirs(os.path.join(options.out_name, chan.name))

        for var in variables:
            h_QCDTemplateFit = f_QCDTemplateFit.Get(f"{sample_QCDTemplateFit.name}_{chan.name}_{var.name}").Clone(
                f"{sample_QCDTemplateFit.name}_{chan.name}_{var.name}_QCD")
            h_MatrixMethod = f_MatrixMethod.Get(f"{sample_MatrixMethod.name}_{chan.name}_{var.name}").Clone(
                f"{sample_MatrixMethod.name}_{chan.name}_{var.name}_MM")
            logging.info(f"Got histograms {h_QCDTemplateFit} {h_MatrixMethod}")

            # Color
            h_QCDTemplateFit.SetLineColor(ROOT.kRed)
            h_QCDTemplateFit.SetMarkerColor(ROOT.kRed)

            # marker size
            h_QCDTemplateFit.SetMarkerSize(0.8)
            h_MatrixMethod.SetMarkerSize(0.8)

            # Rebin
            ratio_range = [0.01, 1.99]
            if "OS-SS" in chan.name:
                h_QCDTemplateFit.Rebin(10)
                h_MatrixMethod.Rebin(10)
                ratio_range = [-3.99, 3.99]

            # check if last plot
            last_plot = var == variables[-1]

            # mc map
            mc_map = {
                sample_QCDTemplateFit: h_QCDTemplateFit,
                sample_MatrixMethod: h_MatrixMethod,
            }

            # canvas
            canv = utils.make_canvas_mc_ratio(mc_map[samples[0]], var, chan, ratio_title=options.ratio_title, x=800, y=800, ratio_range=ratio_range)

            # configure histograms
            canv.configure_histograms(mc_map, options.normalize)

            # top pad
            canv.pad1.cd()
            for s in samples:
                mc_map[s].Draw("pe same")

            # make legend
            canv.make_legend(mc_map, samples)

            # set maximum after creating legend
            canv.set_maximum([mc_map[s] for s in samples], var, mc_map[samples[0]])

            # minimum
            if h_QCDTemplateFit.GetMinimum() < 0 or h_MatrixMethod.GetMinimum() < 0:
                canv.proxy_up.SetMinimum(min(h_QCDTemplateFit.GetMinimum(), h_MatrixMethod.GetMinimum()))

            # bottom pad
            canv.pad2.cd()

            # ratio histograms
            ratios = []
            for i in range(0, len(samples), 2):
                h = mc_map[samples[i]].Clone(f"{mc_map[samples[i]].GetName()}_ratio")
                h.Divide(mc_map[samples[i + 1]])
                ratios += [h]
                h.Draw("same pe")

            # Print out
            canv.print_all(options.analysis_config, chan.name, var.name, multipage_pdf=True,
                           first_plot=first_plot, last_plot=last_plot, as_png=True, logy=False)
            first_plot = False


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-a', '--analysis-config',
                      action="store", dest="analysis_config",
                      help="analysis config file")
    parser.add_option('-n', '--normalize',
                      action="store_true", dest="normalize",
                      help="normalize to luminosity")
    parser.add_option('-t', '--ratio-title',
                      action="store", dest="ratio_title",
                      default="Ratio",
                      help="title of the ratio")
    parser.add_option('--suffix',
                      action="store", dest="suffix",
                      help="suffix for the output name")
    parser.add_option('--stage-out',
                      action="store_true", dest="stage_out",
                      help="copy plots to the www folder")

    # parse input arguments
    options, args = parser.parse_args()

    # analysis configs
    config = options.analysis_config

    # output name
    out_name = config
    if options.suffix:
        out_name = out_name.split("/")
        out_name[0] += "_" + options.suffix
        out_name = "/".join(out_name)
    setattr(options, "out_name", out_name)

    # make output folder if not exist
    if not os.path.isdir(out_name):
        os.makedirs(out_name)

    # do the plotting
    main(options)

    # stage-out to the www folder
    if options.stage_out:
        www.stage_out_plots(out_name, variables, x=300, y=300)
