#!/usr/bin/env python
from charmplot.common import utils
from charmplot.control import channel
from charmplot.control import sample
from charmplot.control import variable
import array
import logging
import os
import re
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

# channel name
channel_names = {
    "mu_minus_Dplus": "W(#mu^{-}#nu)+D^{+}",
    "mu_plus_Dplus": "W(#mu^{+}#nu)+D^{-}",
    "mu_minus_Dstar": "W(#mu^{-}#nu)+D*^{+}",
    "mu_plus_Dstar": "W(#mu^{+}#nu)+D*^{-}",
    "mu_minus_Ds": "W(#mu^{-}#nu)+Ds^{+}",
    "mu_plus_Ds": "W(#mu^{+}#nu)+Ds^{-}",
}

# samples_dict = {
#     "mu_plus_Dplus_LOMG": sample.Sample('mu_plus_Dplus_LOMG', None, **{'add': [], 'subtract': [], 'legendLabel': 'LO MG W+', 'lineColor': 'ROOT.kBlue'}),
#     "mu_minus_Dplus_LOMG": sample.Sample('mu_minus_Dplus_LOMG', None, **{'add': [], 'subtract': [], 'legendLabel': 'LO MG W-', 'lineColor': 'ROOT.kBlue'}),
#     "mu_plus_Dplus_LOMG_pass_w": sample.Sample('mu_plus_Dplus_LOMG_pass_w', None, **{'add': [], 'subtract': [], 'legendLabel': 'LO MG W+', 'lineColor': 'ROOT.kBlue'}),
#     "mu_minus_Dplus_LOMG_pass_w": sample.Sample('mu_minus_Dplus_LOMG_pass_w', None, **{'add': [], 'subtract': [], 'legendLabel': 'LO MG W-', 'lineColor': 'ROOT.kBlue'}),
#     "mu_plus_Dplus_MGFxFx": sample.Sample('mu_plus_Dplus_MGFxFx', None, **{'add': [], 'subtract': [], 'legendLabel': 'MG FxFx W+', 'lineColor': 'ROOT.kRed'}),
#     "mu_minus_Dplus_MGFxFx": sample.Sample('mu_minus_Dplus_MGFxFx', None, **{'add': [], 'subtract': [], 'legendLabel': 'MG FxFx W-', 'lineColor': 'ROOT.kRed'}),
#     "mu_plus_Dplus_MGFxFx_pass_w": sample.Sample('mu_plus_Dplus_MGFxFx_pass_w', None, **{'add': [], 'subtract': [], 'legendLabel': 'MG FxFx W+', 'lineColor': 'ROOT.kRed'}),
#     "mu_minus_Dplus_MGFxFx_pass_w": sample.Sample('mu_minus_Dplus_MGFxFx_pass_w', None, **{'add': [], 'subtract': [], 'legendLabel': 'MG FxFx W-', 'lineColor': 'ROOT.kRed'}),
#     "mu_plus_Dplus_Sherpa2210": sample.Sample('mu_plus_Dplus_Sherpa2210', None, **{'add': [], 'subtract': [], 'legendLabel': 'Sherpa2.2.10 W+', 'lineColor': 'ROOT.kGreen+2'}),
#     "mu_minus_Dplus_Sherpa2210": sample.Sample('mu_minus_Dplus_Sherpa2210', None, **{'add': [], 'subtract': [], 'legendLabel': 'Sherpa2.2.10 W-', 'lineColor': 'ROOT.kGreen+2'}),
#     "mu_plus_Dplus_Sherpa2210_pass_w": sample.Sample('mu_plus_Dplus_Sherpa2210_pass_w', None, **{'add': [], 'subtract': [], 'legendLabel': 'Sherpa2.2.10 W+', 'lineColor': 'ROOT.kGreen+2'}),
#     "mu_minus_Dplus_Sherpa2210_pass_w": sample.Sample('mu_minus_Dplus_Sherpa2210_pass_w', None, **{'add': [], 'subtract': [], 'legendLabel': 'Sherpa2.2.10 W-', 'lineColor': 'ROOT.kGreen+2'}),
# }


# proxy samples
samples_dict = {}
for sampl, color in zip(["LOMG", "MGFxFx", "Sherpa2210"], ["ROOT.kBlue", "ROOT.kRed", "ROOT.kGreen+2"]):
    for channl in channel_names:
        for pass_w in [True, False]:
            name = f"{channl}_{sampl}"
            if pass_w:
                name += "_pass_w"
            samples_dict[name] = sample.Sample(name, None, **{'add': [], 'subtract': [], 'legendLabel': sampl, 'lineColor': color})

print(samples_dict)

# variables to plot
variables = {
    # "D_eta": variable.Variable("D_eta", **{
    #     "label": "#eta(D^{#pm})",
    #     "x_range": [-3.0, 3.0],
    #     "ratio_range": [0.89, 1.11],
    #     "rebin": 2,
    # }),
    "D_pt": variable.Variable("D_pt", **{
        "label": "p_{T}(D^{#pm})",
        "unit": "GeV",
        "x_range": [0, 80],
        "ratio_range": [0.89, 1.11],
        "rebin": 5,
    }),
    # "lep_eta": variable.Variable("lep_eta", **{
    #     "label": "#eta(lep)",
    #     "x_range": [-3.0, 3.0],
    #     "ratio_range": [0.89, 1.11],
    #     "rebin": 2,
    # }),
    # "lep_pt": variable.Variable("lep_pt", **{
    #     "label": "p_{T}(lep)",
    #     "unit": "GeV",
    #     "x_range": [0, 200],
    #     "ratio_range": [0.89, 1.11],
    #     "rebin": 2,
    # }),
    # "met_mt": variable.Variable("met_mt", **{
    #     "label": "m_{T}",
    #     "unit": "GeV",
    #     "x_range": [0, 200],
    #     "ratio_range": [0.89, 1.11],
    #     "rebin": 2,
    # }),
    # "met_met": variable.Variable("met_met", **{
    #     "label": "MET",
    #     "unit": "GeV",
    #     "x_range": [0, 200],
    #     "ratio_range": [0.89, 1.11],
    #     "rebin": 2,
    # }),
}

def add_histograms(histograms, c, var, pass_w = True):
    if pass_w:
        file_name = os.path.join(options.input, f"{c}_{var.name}_pass_w.dat")
    else:
        file_name = os.path.join(options.input, f"{c}_{var.name}.dat")
    f = open(file_name, "r")
    logging.info(f"Succesfully opened file {file_name} as {f}")
    histograms.update(process_file(f, c, var, pass_w))
    f.close()


def create_histogram(variation, x_low, x_high, y, y_up, y_dn, name):
    logging.info(f"creating histogram for {variation}")
    x_bins = x_low + [x_high[-1]]
    x_bins_array = array.array('d', x_bins)
    h = ROOT.TH1D(f"{name}_{variation}", f"{name}_{variation}", len(x_bins) - 1, x_bins_array)
    logging.info(f"created histogram {h}, filling entries now")
    integral = 0
    one_up = 0
    one_dn = 0
    for i in range(h.GetNbinsX()):
        bin_width = h.GetBinWidth(i + 1)
        integral += bin_width * y[i]
        one_up += bin_width * y_up[i]
        one_dn += bin_width * y_dn[i]
        err = (y_up[i] + y_dn[i]) / 2.
        h.SetBinContent(i + 1, y[i])
        h.SetBinError(i + 1, err)
    utils.eprint(f"{variation} {name}: {integral:.3f} + {one_up:.3f} - {one_dn:.3f}")
    return h


def process_file(f, c, var, pass_w):
    histograms = {}
    at_histogram = False
    save_histo = False
    for line in f:
        if "BEGIN HISTO1D" in line:  # noqa: E501
            logging.info(f"{line}")
            variation = re.findall("BEGIN HISTO1D /(.+)", line)  # noqa: W605
            at_histogram = True
            variation = variation[0]
            x_low = []
            x_high = []
            y = []
            y_up = []
            y_dn = []
            print(f"found variation {variation}")
        elif at_histogram and "xlow" in line and "xhigh" in line and "val" in line:
            save_histo = True
        elif "END HISTO1D" in line:
            if at_histogram:
                if not ("[" in variation and "]" in variation):
                    h = create_histogram(variation, x_low, x_high, y, y_up, y_dn, f"{c}_{var.name}")
                    utils.rebin_histogram(h, var)
                    histograms[f"{c}_{variation}{'_pass_w' if pass_w else ''}"] = h
            at_histogram = False
            save_histo = False
        elif save_histo:
            vals = [float(x) for x in line.split()]
            x_low += [vals[0]]
            x_high += [vals[1]]
            y += [vals[2]]
            y_up += [vals[3]]
            y_dn += [vals[4]]
    return histograms


def main(options):

    # save to file
    f = ROOT.TFile(os.path.join(options.output, "output.root"), "RECREATE")
    f.Close()

    for normalize in [True, False]:

        for pass_w in [True, False]:

            for c in options.channels.split(";"):

                # make output folder if not exist
                if not os.path.isdir(os.path.join(options.output, c.split(":")[0])):
                    os.makedirs(os.path.join(options.output, c.split(":")[0]))

                # bookkeeping for later
                first_plot = True

                # vars
                vars_list = options.vars.split(",") if options.vars else [key for key in variables.keys()]
                for v in vars_list:

                    # histograms
                    histograms = {}

                    # var
                    var = variables[v]

                    # check if last plot
                    last_plot = v == vars_list[-1]

                    # check if sub-channels
                    if ":" in c:
                        subchannels = c.split(":")[1].split(",")
                        files = {}
                        for subchannel in subchannels:
                            add_histograms(histograms, subchannel, var, pass_w)
                    else:
                        add_histograms(histograms, c, var, pass_w)

                    # normalize to unity if requested
                    if normalize:
                        for h in histograms.values():
                            h.Scale(1. / h.GetSum())

                    # channel
                    channel_name = c.split(":")[0]
                    chan = channel.Channel(channel_name, [channel_names[channel_name]], "", [], [])

                    # samples
                    samples = [samples_dict[x] for x in histograms.keys() if x in samples_dict]

                    # mc map
                    mc_map = {samples_dict[x]: histograms[x] for x in histograms if x in samples_dict}

                    # canvas
                    yaxis_label = f"d#sigma / d{var.label} [pb]"
                    if normalize:
                        yaxis_label = "Normalized Entries"
                    canv = utils.make_canvas_mc_ratio(mc_map[samples[0]], var, chan, "Ratio", x=800, y=800, events=yaxis_label)

                    # configure histograms
                    canv.configure_histograms(mc_map, True)

                    # top pad
                    canv.pad1.cd()
                    errors = []
                    canv.pad1.cd()
                    for s in samples:
                        fcolor = mc_map[s].GetLineColor()
                        gr_mc_stat_err, _ = utils.make_stat_err(mc_map[s])
                        gr_mc_stat_err.SetLineColor(fcolor)
                        gr_mc_stat_err.SetFillColorAlpha(fcolor, 0.25)
                        gr_mc_stat_err.SetFillStyle(1001)
                        errors += [gr_mc_stat_err]
                        gr_mc_stat_err.Draw("e2")
                        mc_map[s].Draw("hist same")

                    # make legend
                    canv.make_legend(mc_map, samples, print_yields=(not normalize), show_error=False, yields_unit="pb", leg_offset=-0.05)

                    # set maximum after creating legend
                    canv.set_maximum([mc_map[s] for s in samples], var, mc_map[samples[0]])

                    # bottom pad
                    canv.pad2.cd()

                    # ratio histograms
                    ratios = []
                    denominator = mc_map[samples[0]].Clone(f"{mc_map[samples[0]].GetName()}_denominator")
                    for i in range(0, denominator.GetNbinsX() + 2):
                        denominator.SetBinError(i, 0)
                    for i in range(0, len(samples)):
                        h = mc_map[samples[i]].Clone(f"{mc_map[samples[i]].GetName()}_ratio")
                        h.Divide(denominator)
                        ratios += [h]
                        fcolor = mc_map[samples[i]].GetLineColor()
                        gr_mc_stat_err, _ = utils.make_stat_err(h)
                        gr_mc_stat_err.SetLineColor(fcolor)
                        gr_mc_stat_err.SetFillColorAlpha(fcolor, 0.25)
                        gr_mc_stat_err.SetFillStyle(1001)
                        errors += [gr_mc_stat_err]
                        gr_mc_stat_err.Draw("e2")
                        h.Draw("hist same")
                    
                    # ratio range
                    if normalize:
                        canv.proxy_dn.SetMaximum(1.19)
                        canv.proxy_dn.SetMinimum(0.81)

                    # Print out
                    canv.print_all(options.output, channel_name + ("_pass_w" if pass_w else "") + ("_norm" if normalize else ""),
                                   v, multipage_pdf=True, first_plot=first_plot, last_plot=last_plot, as_png=False)
                    first_plot = False


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-i', '--input',
                      action="store", dest="input",
                      help="analysis config file")
    parser.add_option('-o', '--output',
                      action="store", dest="output",
                      help="output folder name",
                      default="rivet_plots")
    parser.add_option('-c', '--channels',
                      action="store", dest="channels",
                      default="mu_minus_Dplus;mu_plus_Dplus;mu_minus_Dstar;mu_plus_Dstar;mu_minus_Ds;mu_plus_Ds",
                      help="run over a subset of channels")
    parser.add_option('-v', '--vars',
                      action="store", dest="vars",
                      help="run over a subset of variables (comma separated)")

    # parse input arguments
    options, args = parser.parse_args()

    # output name
    out_name = options.output

    # make output folder if not exist
    if not os.path.isdir(out_name):
        os.makedirs(out_name)

    # do the plotting
    main(options)
