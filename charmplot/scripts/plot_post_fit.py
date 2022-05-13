#!/usr/bin/env python
from charmplot.common import utils
from charmplot.common import www
from charmplot.control.channel import Channel
from charmplot.control import globalConfig
from charmplot.control import tools
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


def modify_sample_name(sample_name, channel):
    # if "WplusD_Matched" in sample_name:
    #     if "minus" in channel.name:
    #         sample_name += "_minus"
    #     elif "plus" in channel.name:
    #         sample_name += "_plus"
    return sample_name


def get_err_hist(f, par, variation, default):
    h_err = None
    name = default.split("postFit")[0]
    h_temp = f.Get(f"{name}{par}_{variation}_postFit")
    if h_temp:
        h_err = h_temp.Clone(f"{h_temp.GetName()}_err_{variation}")
    else:
        # logging.warning(f"Histogram {name}{par}_{variation}_postFit not found in file {os.path.basename(f.GetName())}. Using {default} instead.")
        h_mc_default = f.Get(default)
        if h_mc_default:
            h_err = h_mc_default.Clone(f"{h_mc_default.GetName()}_err_{variation}")
    return h_err


def main(options, conf):
    trex_histogram_folder = os.path.join(options.trex_input, "Histograms")

    # out root file
    outfile = ROOT.TFile(f"post_fit/{conf.out_name}/outpoot.root", "RECREATE")

    # channels
    channels = []

    # fitted variable
    var = conf.get_var(options.var)
    logging.info(f"Got variable {var}")

    # parse input
    for file in os.listdir(trex_histogram_folder):
        if file.endswith("postFit.root"):
            channel_name = file.replace("_postFit.root", "")
            logging.info(f"Searching for channel {channel_name}...")
            channel = conf.get_channel(channel_name)
            if not channel:
                logging.warning(f"Channel not found for string {channel_name}")
            else:
                logging.info(f"Found channel {channel_name}")
                if not options.skip_channel in channel.name:
                    channels += [channel]

    # sort channels
    individual_plots = []
    OS_minus_SS_plots = []
    OS_minus_SS_total = {'+': [], '-': []}
    OS_minus_SS_total_minus = {'+': [], '-': []}
    OS_minus_SS_total_plus = {'+': [], '-': []}
    for channel in channels:
        individual_plots += [{'+': [channel], '-': []}]
        logging.info(f"Added channel {channel.name}..")
    for channel_OS in channels:
        if "OS_" not in channel_OS.name:
            continue
        for channel_SS in channels:
            if channel_SS.name == channel_OS.name.replace("OS_", "SS_"):
                OS_minus_SS_plots += [{'+': [channel_OS], '-': [channel_SS]}]
                if "0tag" in channel_SS.name:
                    OS_minus_SS_total['+'] += [channel_OS]
                    OS_minus_SS_total['-'] += [channel_SS]
                    if "minus" in channel_SS.name:
                        OS_minus_SS_total_minus['+'] += [channel_OS]
                        OS_minus_SS_total_minus['-'] += [channel_SS]
                    elif "plus" in channel_SS.name:
                        OS_minus_SS_total_plus['+'] += [channel_OS]
                        OS_minus_SS_total_plus['-'] += [channel_SS]
                break

    # get correlation matrix
    logging.info("Loading correlation matrix...")
    corr_dict = tools.parse_yaml_file(os.path.join(options.trex_input, "CorrelationMatrix.yaml"))
    corr_parameters = None
    corr_correlation_rows = None
    for x in corr_dict:
        if 'parameters' in x:
            corr_parameters = x['parameters']
        elif 'correlation_rows' in x:
            corr_correlation_rows = x['correlation_rows']
    n_pars = len(corr_parameters)

    # plots = individual_plots + OS_minus_SS_plots + [OS_minus_SS_total]
    # plots = individual_plots + OS_minus_SS_plots
    plots = OS_minus_SS_plots + [OS_minus_SS_total, OS_minus_SS_total_minus, OS_minus_SS_total_plus]
    # plots = [OS_minus_SS_total]

    for plot in plots:

        # Skip channel if list in dict empty
        if not plot['+']:
            continue

        # create channel
        channel_temp = plot['+'][0]
        channels_all = plot['+'] + plot['-']
        channel_name = "_".join([channel.name for channel in channels_all])

        # inclusive?
        minus = False
        plus = False
        inclusive = False
        inclusive_minus = False
        inclusive_plus = False
        if len(channels_all) > 2:
            for c in channels_all:
                if "minus" in c.name:
                    minus = True
                elif "plus" in c.name:
                    plus = True
            if minus and plus:
                inclusive = True
            elif minus and not plus:
                inclusive_minus = True
            elif plus and not minus:
                inclusive_plus = True

        # labels
        labels = []
        for label in channel_temp.label:
            if inclusive and "channel" in label:
                label = "W^{#pm} channel"
            if len(channels_all) > 2 and "bin" in label:
                label = "inclusive 0-tag"
            if len(plot['-']) > 0:
                label = label.replace("OS", "OS-SS")
            if len(plot['+']) > 1 and "tag" in label:
                label = label.split(",")[0]
            labels += [label]
        labels[-1] += ', post-fit'

        if inclusive:
            channel_name = "0tag_inclusive"
        elif inclusive_minus:
            channel_name = "0tag_inclusive_minus"
        elif inclusive_plus:
            channel_name = "0tag_inclusive_plus"

        chan = Channel(channel_name, labels, channel_temp.lumi, [], [])

        # read files
        file_suffix = ""
        if var.name not in ["Dmeson_m", "Dmeson_mdiff", "Dmeson_m_fit", "Dmeson_mdiff_fit"]:
            file_suffix = f"{var.name.replace('Dmeson', '')}"
        files = {}
        for channel in channels_all:
            files[channel] = ROOT.TFile(os.path.join(trex_histogram_folder, f"{channel.name}{file_suffix}_postFit.root"))

        # samples
        sample_names = []
        samples = []
        for channel in channels_all:
            for sample_name in channel.samples:
                sample = conf.get_sample(sample_name)
                if sample.shortName not in sample_names:
                    sample_names.append(sample.shortName)
                    samples.append(sample)

        # get mc samples
        h_bkg = None
        mc_map = {}
        for sample in samples:
            h_sum = None
            sample_name = sample.shortName

            # positive channels
            for channel in plot['+']:

                h_temp = files[channel].Get(f"h_{modify_sample_name(sample_name, channel)}_postFit")
                # print(files[channel].GetName()," ",h_temp," ",modify_sample_name(sample_name, channel)," ",channel.name)
                if h_temp:
                    if h_sum is None:
                        h_sum = h_temp.Clone(f"{h_temp.GetName()}_{chan.name}")
                    else:
                        h_sum.Add(h_temp)
                    if options.subtract_background and options.signal not in modify_sample_name(sample_name, channel):
                        if h_bkg is None:
                            h_bkg = h_temp.Clone(f"{chan.name}_bkg")
                        else:
                            h_bkg.Add(h_temp)

            # negative channels
            for channel in plot['-']:
                h_temp = files[channel].Get(f"h_{modify_sample_name(sample_name, channel)}_postFit")
                if h_temp:
                    if h_sum is None:
                        h_sum = h_temp.Clone(f"{h_temp.GetName()}_{chan.name}")
                        h_sum.Scale(-1.)
                    else:
                        h_sum.Add(h_temp, -1)
                    if options.subtract_background and options.signal not in modify_sample_name(sample_name, channel):
                        if h_bkg is None:
                            h_bkg = h_temp.Clone(f"{chan.name}_bkg")
                            h_bkg.Scale(-1.)
                        else:
                            h_bkg.Add(h_temp, -1)

            if h_sum and abs(h_sum.GetSum()) > 1e-2:
                if (options.subtract_background and options.signal in modify_sample_name(sample_name, channel)) or not options.subtract_background:
                    mc_map[sample] = h_sum

        # get data
        h_data = None
        for channel in plot['+']:
            h_temp = files[channel].Get("h_Data")
            if not h_data:
                h_data = h_temp.Clone(f"{h_temp.GetName()}_{chan.name}")
            else:
                h_data.Add(h_temp)
            if len(plot['-']):
                for channel_SS in channels:
                    if channel_SS.name == channel.name.replace("OS_", "SS_"):
                        h_temp_SS = files[channel_SS].Get("h_Data")
                        h_data.Add(h_temp_SS, -1)
        if options.subtract_background and h_bkg:
            h_data.Add(h_bkg, -1)

        # get mc tot
        h_mc_tot = None
        for channel in plot['+']:
            h_temp = files[channel].Get("h_tot_postFit")
            if not h_mc_tot:
                h_mc_tot = h_temp.Clone(f"{h_temp.GetName()}_{chan.name}")
            else:
                h_mc_tot.Add(h_temp)
            if len(plot['-']):
                for channel_SS in channels:
                    if channel_SS.name == channel.name.replace("OS_", "SS_"):
                        h_temp_SS = files[channel_SS].Get("h_tot_postFit")
                        h_mc_tot.Add(h_temp_SS, -1)
        if options.subtract_background and h_bkg:
            h_mc_tot.Add(h_bkg, -1)

        # canvas
        canv = utils.make_canvas(h_data, var, chan, x=800, y=800)

        # configure histograms
        canv.configure_histograms(mc_map, h_data, style=conf.style)

        # stack and total mc
        hs = utils.make_stack(samples, mc_map)

        # stat error
        g_mc_tot_err, g_mc_tot_err_only = utils.make_stat_err(h_mc_tot)

        # systematic uncertainties
        h_mc_tot_err_histograms_Up = []
        h_mc_tot_err_histograms_Dn = []
        for par in corr_parameters:
            h_mc_tot_Up = None
            h_mc_tot_Dn = None
            for channel in plot['+']:
                h_temp_up = get_err_hist(files[channel], par, "Up", "h_tot_postFit")
                h_temp_dn = get_err_hist(files[channel], par, "Down", "h_tot_postFit")
                if h_temp_up:
                    if not h_mc_tot_Up:
                        h_mc_tot_Up = h_temp_up.Clone(f"{h_temp_up.GetName()}_{chan.name}_err_up")
                    else:
                        h_mc_tot_Up.Add(h_temp_up)
                if h_temp_dn:
                    if not h_mc_tot_Dn:
                        h_mc_tot_Dn = h_temp_dn.Clone(f"{h_temp_dn.GetName()}_{chan.name}_err_dn")
                    else:
                        h_mc_tot_Dn.Add(h_temp_dn)
            for channel in plot['-']:
                h_temp_up = get_err_hist(files[channel], par, "Up", "h_tot_postFit")
                h_temp_dn = get_err_hist(files[channel], par, "Down", "h_tot_postFit")
                if h_temp_up:
                    if not h_mc_tot_Up:
                        h_mc_tot_Up = h_temp_up.Clone(f"{h_temp_up.GetName()}_{chan.name}_err_up")
                        h_mc_tot_Up.Scale(-1.)
                    else:
                        h_mc_tot_Up.Add(h_temp_up, -1)
                if h_temp_dn:
                    if not h_mc_tot_Dn:
                        h_mc_tot_Dn = h_temp_dn.Clone(f"{h_temp_dn.GetName()}_{chan.name}_err_dn")
                        h_mc_tot_Dn.Scale(-1)
                    else:
                        h_mc_tot_Dn.Add(h_temp_dn, -1)

            # subtract nominal
            if options.subtract_background and h_bkg:
                h_mc_tot_Up.Add(h_bkg, -1)
                h_mc_tot_Dn.Add(h_bkg, -1)
            h_mc_tot_Up.Add(h_mc_tot, -1)
            h_mc_tot_Dn.Add(h_mc_tot, -1)
            h_mc_tot_err_histograms_Up += [h_mc_tot_Up]
            h_mc_tot_err_histograms_Dn += [h_mc_tot_Dn]

        # calculate error bands
        for x in range(1, h_mc_tot.GetNbinsX() + 1):
            # add stat uncertaitny from the histogram
            if file_suffix:
                g_mc_tot_err.GetEYlow()[x - 1] = g_mc_tot_err.GetEYlow()[x - 1]**2
                g_mc_tot_err.GetEYhigh()[x - 1] = g_mc_tot_err.GetEYhigh()[x - 1]**2
            else:
                g_mc_tot_err.GetEYlow()[x - 1] = 0.
                g_mc_tot_err.GetEYhigh()[x - 1] = 0.
            # off diagonal
            for i in range(n_pars):
                for j in range(i):
                    if h_mc_tot_err_histograms_Up[i] and h_mc_tot_err_histograms_Up[j] and h_mc_tot_err_histograms_Dn[i] and h_mc_tot_err_histograms_Dn[j]:
                        corr = corr_correlation_rows[i][j]
                        err_up_i = h_mc_tot_err_histograms_Up[i].GetBinContent(x)
                        err_dn_i = h_mc_tot_err_histograms_Dn[i].GetBinContent(x)
                        err_up_j = h_mc_tot_err_histograms_Up[j].GetBinContent(x)
                        err_dn_j = h_mc_tot_err_histograms_Dn[j].GetBinContent(x)
                        err_up = err_up_i * err_up_j * corr * 2
                        err_dn = err_dn_i * err_dn_j * corr * 2
                        g_mc_tot_err.GetEYhigh()[x - 1] += err_up
                        g_mc_tot_err.GetEYlow()[x - 1] += err_dn
            # diagonal
            for i in range(n_pars):
                if h_mc_tot_err_histograms_Up[i] and h_mc_tot_err_histograms_Dn[i]:
                    err_up_i = h_mc_tot_err_histograms_Up[i].GetBinContent(x)
                    err_dn_i = h_mc_tot_err_histograms_Dn[i].GetBinContent(x)
                    err_up = err_up_i * err_up_i
                    err_dn = err_dn_i * err_dn_i
                    g_mc_tot_err.GetEYhigh()[x - 1] += err_up
                    g_mc_tot_err.GetEYlow()[x - 1] += err_dn

            # final
            g_mc_tot_err.GetEYhigh()[x - 1] = math.sqrt(g_mc_tot_err.GetEYhigh()[x - 1])
            g_mc_tot_err.GetEYlow()[x - 1] = math.sqrt(g_mc_tot_err.GetEYlow()[x - 1])
            if h_mc_tot.GetBinContent(x) > 0:
                g_mc_tot_err_only.GetEYhigh()[x - 1] = g_mc_tot_err.GetEYhigh()[x - 1] / h_mc_tot.GetBinContent(x)
                g_mc_tot_err_only.GetEYlow()[x - 1] = g_mc_tot_err.GetEYlow()[x - 1] / h_mc_tot.GetBinContent(x)
            else:
                g_mc_tot_err_only.GetEYhigh()[x - 1] = 0
                g_mc_tot_err_only.GetEYlow()[x - 1] = 0

        # ratio
        h_ratio = utils.make_ratio(h_data, h_mc_tot)

        # data graph
        gr_data = utils.get_gr_from_hist(h_data)
        gr_ratio = utils.get_gr_from_hist(h_ratio)

        # top pad
        canv.pad1.cd()

        # make legend
        canv.make_legend(h_data, h_mc_tot, mc_map, samples, print_yields=True, show_error=False)

        # normalize bins to unity
        if var.per_unit:
            utils.normalize_to_unit(hs, hists=[h_mc_tot, h_data], grs=[g_mc_tot_err, gr_data])

        hs.Draw("same hist")
        h_mc_tot.Draw("same hist")
        g_mc_tot_err.Draw("e2")
        gr_data.Draw("pe0")

        # set maximum after creating legend
        canv.set_maximum((h_data, h_mc_tot), var, mc_min=utils.get_mc_min(mc_map, samples))

        # find minimum
        if '0tag' in channel.name:
            min_negative = {}
            for s in samples:
                if s not in mc_map:
                    continue
                h = mc_map[s]
                for i in range(1, h.GetNbinsX() + 1):
                    if h.GetBinContent(i) < 0:
                        if i not in min_negative:
                            min_negative[i] = 0
                        min_negative[i] += h.GetBinContent(i)
            if min_negative.values():
                min_negative = min(min_negative.values())
                canv.proxy_up.SetMinimum(min_negative)

        # bottom pad
        ROOT.gPad.RedrawAxis()
        canv.pad2.cd()
        g_mc_tot_err_only.Draw("le2")
        gr_ratio.Draw("pe0")

        # ratio range
        ratio_range = [1.01 - float(options.ratio_range) / 100., 0.99 + float(options.ratio_range) / 100.]
        canv.set_ratio_range(ratio_range[0], ratio_range[1], override=True)
        ROOT.gPad.RedrawAxis()

        # Print out
        canv.print(f"post_fit/{conf.out_name}/{chan.name}_{var.name}.pdf")
        outfile.cd()
        gr_data.Write(f"{h_mc_tot.GetName()}_data")
        h_mc_tot.Write()

        logging.info(f"finished processing channel {channel.name}")

    # close output file
    outfile.Close()


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-a', '--analysis-config',
                      action="store", dest="analysis_config",
                      help="analysis config file")
    parser.add_option('-b', '--subtract-background',
                      action="store_true", dest="subtract_background",
                      help="subtract background from data and total MC")
    parser.add_option('-s', '--signal',
                      action="store", dest="signal",
                      default="WplusD",
                      help="name of the signal sample (to be used with '-b')")
    parser.add_option('-v', '--var',
                      action="store", dest="var",
                      help="fitted variable",
                      default="Dmeson_m")
    parser.add_option('-r', '--ratio-range',
                      action="store", dest="ratio_range",
                      default=50,
                      help="range of the ratio y-axis")
    parser.add_option('-k', '--skip-channel',
                      action="store", dest="skip_channel",
                      default="1tag",
                      help="skip channels that include this in the name")
    parser.add_option('--suffix',
                      action="store", dest="suffix",
                      help="suffix for the output name")
    parser.add_option('--stage-out',
                      action="store_true", dest="stage_out",
                      help="copy plots to the www folder")
    parser.add_option('--trex-input',
                      action="store", dest="trex_input",
                      help="import post-fit trex plots")
    parser.add_option('-t', '--threads',
                      action="store", dest="threads",
                      help="number of threads",
                      default=8)

    # parse input arguments
    options, args = parser.parse_args()

    # analysis configs
    config = options.analysis_config

    # output name
    out_name = config.replace(".yaml", "").replace(".yml", "")
    if options.suffix:
        out_name = out_name.split("/")
        if out_name[0] != ".":
            out_name[0] += "_" + options.suffix
            out_name = "/".join(out_name)
        else:
            out_name[1] += "_" + options.suffix
            out_name = "/".join(out_name)

    # config object
    conf = globalConfig.GlobalConfig(config, out_name)

    # make output folder if not exist
    if not os.path.isdir(os.path.join("post_fit", out_name)):
        os.makedirs(os.path.join("post_fit", out_name))

    # do the plotting
    main(options, conf)

    # stage-out to the www folder
    if options.stage_out:
        www.stage_out_plots(conf.out_name, conf.get_variables(), x=300, y=300)
