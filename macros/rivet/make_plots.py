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

# samples
SAMPLES = ["LOMG", "NLOMG[_VarBand]", "QCDScale[_VarBand]", "MGFxFx", "MGFxFxReshower", "Sherpa2210", "Sherpa2211"]
COLORS = ["ROOT.kBlue", "ROOT.kBlack", "ROOT.kRed", "ROOT.kYellow+2", "ROOT.kBlack", "ROOT.kGreen+2", "ROOT.kViolet"]

# channel name
channel_names = {
    "mu_minus_Dplus": "W(#mu^{-}#nu)+D^{+}",
    "mu_plus_Dplus": "W(#mu^{+}#nu)+D^{-}",
    "mu_minus_Dstar": "W(#mu^{-}#nu)+D*^{+}",
    "mu_plus_Dstar": "W(#mu^{+}#nu)+D*^{-}",
    "mu_minus_Ds": "W(#mu^{-}#nu)+Ds^{+}",
    "mu_plus_Ds": "W(#mu^{+}#nu)+Ds^{-}",
    "mu_minus_LambdaC": "W(#mu^{-}#nu)+#Lambda_{C}^{+}",
    "mu_plus_LambdaC": "W(#mu^{+}#nu)+#Lambda_{C}^{-}",
    "mu_minus_OmegaC": "W(#mu^{-}#nu)+#Omega_{C}^{+}",
    "mu_plus_OmegaC": "W(#mu^{+}#nu)+#Omega_{C}^{-}",
    "mu_minus_XiCplus": "W(#mu^{-}#nu)+#Xi_{C}^{+}",
    "mu_plus_XiCplus": "W(#mu^{+}#nu)+#Xi_{C}^{-}",
    "mu_minus_XiCzero": "W(#mu^{-}#nu)+#Xi_{C}0",
    "mu_plus_XiCzero": "W(#mu^{+}#nu)+#Xi_{C}0",

    "mu_minus_mesons": "W(#mu^{-}#nu)+c-meson^{+}",
    "mu_plus_mesons": "W(#mu^{+}#nu)+c-meson^{-}",

    "mu_Dplus": "W(#mu#nu)+D",
    "mu_Dstar": "W(#mu#nu)+D*",
    "mu_Dzero": "W(#mu#nu)+D0",
    "mu_Ds": "W(#mu#nu)+Ds",
    "mu_Dmeson": "W(#mu#nu)+Dmeson",
    "mu_Baryon": "W(#mu#nu)+Baryon",
}

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
        "x_range": [8, 98],
        "ratio_range": [0.89, 1.11],
        "bins": [8, 12, 20, 40, 80, 98],
        "logx": True,
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


def add_histograms(histograms, c, var, pass_w=True, normalize=True):
    if pass_w:
        file_name = os.path.join(options.input, f"{c}_{var.name}_pass_w.dat")
    else:
        file_name = os.path.join(options.input, f"{c}_{var.name}.dat")
    f = open(file_name, "r")
    logging.info(f"Succesfully opened file {file_name} as {f}")
    histograms.update(process_file(f, c, var, pass_w, normalize))
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


def process_file(f, c, var, pass_w, normalize):
    histograms = {}
    at_histogram = False
    save_histo = False
    for line in f:
        if "BEGIN HISTO1D" in line:
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
                # if not ("[" in variation and "]" in variation):
                if variation in SAMPLES:
                    h = create_histogram(variation, x_low, x_high, y, y_up, y_dn,
                                         f"{c}_{var.name}{'_pass_w' if pass_w else ''}{'_normalized' if normalize else ''}")
                    utils.rebin_histogram(h, var)
                    if var.bins:
                        h = h.Rebin(len(var.bins) - 1, f"{h.GetName()}_rebin", array.array('d', var.bins))
                        varN = h.GetBinContent(len(var.bins) - 1)
                        varN1 = h.GetBinContent(len(var.bins))
                        errN = h.GetBinError(len(var.bins) - 1)
                        errN1 = h.GetBinError(len(var.bins))
                        h.SetBinContent(len(var.bins) - 1, varN + varN1)
                        h.SetBinError(len(var.bins) - 1, (errN**2 + errN1**2)**(0.5))
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

    # proxy samples
    samples_dict = {}
    for sampl, color in zip(SAMPLES, COLORS):
        if sampl not in options.samples.split(","):
            continue
        for channl in channel_names:
            for pass_w in [True, False]:
                name = f"{channl}_{sampl}"
                if pass_w:
                    name += "_pass_w"
                legendName = sampl.split("[_VarBand]")[0]
                samples_dict[name] = sample.Sample(name, None, **{'add': [], 'subtract': [], 'legendLabel': legendName, 'lineColor': color})

    print(samples_dict)

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
                        for subchannel in subchannels:
                            add_histograms(histograms, subchannel, var, pass_w, normalize)
                    else:
                        add_histograms(histograms, c, var, pass_w, normalize)

                    # channel
                    channel_name = c.split(":")[0]
                    chan = channel.Channel(channel_name, [channel_names[channel_name]], "", [], [])

                    # samples
                    if ":" not in c:
                        samples = [samples_dict[x] for x in histograms.keys() if x in samples_dict]
                        mc_map = {samples_dict[x]: histograms[x] for x in histograms if x in samples_dict}
                    else:
                        samples = []
                        mc_map = {}

                        for sampl in SAMPLES:
                            if sampl not in options.samples.split(","):
                                continue
                            h_sum = None
                            for h in [x for x in histograms if x.replace('_pass_w', '').replace('_normalized', '').endswith(sampl)]:
                                print(h)
                                if not h_sum:
                                    h_sum = histograms[h].Clone(f"{histograms[h].GetName()}_{channel_name}")
                                else:
                                    h_sum.Add(histograms[h])
                            s = samples_dict[f"{channel_name}_{sampl}{'_pass_w' if pass_w else ''}"]
                            samples += [s]
                            if h_sum:
                                mc_map[s] = h_sum

                    # normalize to unity if requested
                    if normalize:
                        for s in mc_map:
                            mc_map[s].Scale(1. / mc_map[s].GetSum())

                    # canvas
                    yaxis_label = f"d#sigma / d{var.label} [pb]"
                    if normalize:
                        yaxis_label = "Normalized Entries"

                    # first available histo
                    for s in samples:
                        if s in mc_map and mc_map[s]:
                            h0 = mc_map[s]
                            break

                    canv = utils.make_canvas_mc_ratio(h0, var, chan, "Ratio", x=800, y=800, events=yaxis_label)

                    # configure histograms
                    canv.configure_histograms(mc_map, True)

                    # logx
                    if var.logx:
                        canv.pad1.SetLogx()
                        canv.pad2.SetLogx()
                        canv.proxy_dn.GetXaxis().SetMoreLogLabels()
                        canv.proxy_dn.GetXaxis().SetLabelOffset(-0.03)

                    # top pad
                    errors = []
                    canv.pad1.cd()
                    for s in samples:
                        fcolor = mc_map[s].GetLineColor()
                        gr_mc_stat_err, _ = utils.make_stat_err(mc_map[s])
                        gr_mc_stat_err.SetLineColor(fcolor)
                        # gr_mc_stat_err.SetFillColorAlpha(fcolor, 0.25)
                        gr_mc_stat_err.SetFillStyle(3345)
                        errors += [gr_mc_stat_err]
                        gr_mc_stat_err.Draw("e2")
                        mc_map[s].Draw("hist same")

                    # make legend
                    canv.make_legend(mc_map, samples, print_yields=(not normalize), show_error=False, yields_unit="pb", leg_offset=-0.05)

                    # set maximum after creating legend
                    canv.set_maximum([mc_map[s] for s in samples], var, mc_map[samples[0]])

                    # bottom pad
                    ROOT.gPad.RedrawAxis()
                    canv.pad2.cd()

                    # output file
                    f = ROOT.TFile(os.path.join(options.output, "output.root"), "UPDATE")
                    f.cd()

                    # ratio histograms
                    ratios = []
                    denominator = mc_map[samples[0]].Clone(f"{mc_map[samples[0]].GetName()}_denominator")
                    for i in range(0, denominator.GetNbinsX() + 2):
                        denominator.SetBinError(i, 0)
                    for i in range(0, len(samples)):
                        h = mc_map[samples[i]].Clone(f"{mc_map[samples[i]].GetName()}_ratio")
                        h.Divide(denominator)
                        if normalize and not pass_w:
                            h.Write(f"{channel_name}_{samples[i].name}")
                        ratios += [h]
                        fcolor = mc_map[samples[i]].GetLineColor()
                        gr_mc_stat_err, _ = utils.make_stat_err(h)
                        gr_mc_stat_err.SetLineColor(fcolor)
                        # gr_mc_stat_err.SetFillColorAlpha(fcolor, 0.25)
                        gr_mc_stat_err.SetFillStyle(3345)
                        errors += [gr_mc_stat_err]
                        gr_mc_stat_err.Draw("e2")
                        h.Draw("hist same")

                    # ratio range
                    if normalize:
                        canv.proxy_dn.SetMaximum(1.19)
                        canv.proxy_dn.SetMinimum(0.81)

                    # Print out
                    ROOT.gPad.RedrawAxis()
                    canv.print_all(options.output, channel_name + ("_pass_w" if pass_w else "") + ("_norm" if normalize else ""),
                                   v, multipage_pdf=True, first_plot=first_plot, last_plot=last_plot, as_png=False)
                    first_plot = False

                    # close file
                    f.Close()


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
                      default="mu_Dmeson:mu_minus_Dplus,mu_plus_Dplus,mu_minus_Dzero,mu_plus_Dzero,mu_minus_Ds,mu_plus_Ds;mu_Dplus:mu_minus_Dplus,mu_plus_Dplus;mu_Dstar:mu_minus_Dstar,mu_plus_Dstar;mu_Dzero:mu_minus_Dzero,mu_plus_Dzero;mu_Ds:mu_minus_Ds,mu_plus_Ds;mu_Baryon:mu_minus_LambdaC,mu_plus_LambdaC,mu_minus_OmegaC,mu_plus_OmegaC,mu_minus_XiCplus,mu_plus_XiCplus,mu_minus_XiCzero,mu_plus_XiCzero",
                      help="run over a subset of channels")
    parser.add_option('-s', '--samples',
                      action="store", dest="samples",
                      default="LOMG,NLOMG[_VarBand],QCDScale[_VarBand],MGFxFx,MGFxFxReshower,Sherpa2210,Sherpa2211",
                      help="the samples to run over")
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
