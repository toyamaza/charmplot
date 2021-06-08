#!/usr/bin/env python
from charmplot.common import utils
from charmplot.common import www
from charmplot.control.channel import Channel
from charmplot.control import globalConfig
from charmplot.control import tools
import logging
import math
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


def get_err_hist(f, par, variation, default):
    h_err = None
    name = default.split("postFit")[0]
    h_temp = f.Get(f"{name}{par}_{variation}_postFit")
    if h_temp:
        h_err = h_temp.Clone(f"{h_temp.GetName()}_err_{variation}")
    else:
        logging.warning(f"Histogram not found for variation {par}_{variation} in file {os.path.basename(f.GetName())}. Using {default} instead.")
        h_mc_default = f.Get(default)
        if h_mc_default:
            h_err = h_mc_default.Clone(f"{h_mc_default.GetName()}_err_{variation}")
    return h_err


def main(options, conf):
    trex_histogram_folder = os.path.join(options.trex_input, "Histograms")

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
                logging.error(f"Channel not found for string {channel_name}")
            else:
                logging.info(f"Found channel {channel_name}")
            channels += [channel]

    # sort channels
    individual_plots = []
    OS_minus_SS_plots = []
    OS_minus_SS_total = {'+': [], '-': []}
    for channel in channels:
        individual_plots += [{'+': [channel], '-': []}]
        logging.info(f"Added channel {channel.name}..")
    for channel_OS in channels:
        if "OS_" not in channel_OS.name:
            continue
        # if "1tag" not in channel_OS.name:
        #     continue
        for channel_SS in channels:
            if channel_SS.name == channel_OS.name.replace("OS_", "SS_"):
                OS_minus_SS_plots += [{'+': [channel_OS], '-': [channel_SS]}]
                if "0tag" in channel_SS.name:
                    OS_minus_SS_total['+'] += [channel_OS]
                    OS_minus_SS_total['-'] += [channel_SS]
                break

    # get correlation matrix
    logging.info(f"Loading correlation matrix...")
    corr_dict = tools.parse_yaml_file(os.path.join(options.trex_input, "CorrelationMatrix.yaml"))
    corr_parameters = None
    corr_correlation_rows = None
    for x in corr_dict:
        if 'parameters' in x:
            corr_parameters = x['parameters']
        elif 'correlation_rows' in x:
            corr_correlation_rows = x['correlation_rows']
    n_pars = len(corr_parameters)

    plots = individual_plots + OS_minus_SS_plots + [OS_minus_SS_total]
    # plots = individual_plots + OS_minus_SS_plots
    # plots = OS_minus_SS_plots

    print(f"individual_plots {individual_plots}")
    print(f"OS_minus_SS_plots {OS_minus_SS_plots}")

    subtract_manually = True

    for plot in plots:

        # create channel
        channel_temp = plot['+'][0]
        channels_all = plot['+'] + plot['-']
        channel_name = "_".join([channel.name for channel in channels_all])
        if len(channels_all) > 2:
            channel_name = "0tag_inclusive"

        # labels
        labels = []
        for label in channel_temp.label:
            if len(plot['-']) > 0:
                label = label.replace("OS", "OS-SS")
            if len(plot['+']) > 1 and "tag" in label:
                label = label.split(",")[0]
            labels += [label]
        labels[-1] += ', post-fit'
        chan = Channel(channel_name, labels, channel_temp.lumi, [], [])

        # read files
        files = {}
        for channel in channels_all:
            files[channel] = ROOT.TFile(os.path.join(trex_histogram_folder, f"{channel.name}_postFit.root"))

        # samples
        sample_names = []
        samples = []
        for channel in channels_all:
            for sample_name in channel.samples:
                sample = conf.get_sample(sample_name)
                if sample.shortName not in sample_names:
                    sample_names.append(sample.shortName)
                    samples.append(sample)

        # put MockMC last
        if len(plot['-']) and 'MockMC' in samples[0].shortName:
            samples += [samples.pop(0)]
            sample_names += [sample_names.pop(0)]

        # get mc samples
        mc_map = {}
        for sample in samples:
            h_sum = None

            # positive channels
            for channel in plot['+']:

                # channel charge
                channel_charge = ""
                if options.plus_minus_ratio and ('MockMC' in sample.shortName or 'Wjets_emu_Matched' in sample.shortName):
                    if "minus" in channel.name:
                        channel_charge = "_minus"
                    elif "plus" in channel.name:
                        channel_charge = "_plus"

                # pt bin
                pt_bin = ""
                if 'MockMC' in sample.shortName or 'Wjets_emu_Matched' in sample.shortName:
                    pt_bin = "_" + re.findall("(pt_bin[0-9])", channel.name)[0]

                # truth pt bin
                truth_pt_bin = ""
                if pt_bin and 'Wjets_emu_Matched' in sample.shortName and "truth_pt" not in sample.shortName:
                    truth_pt_bin = f"_truth{pt_bin}"

                # mockMC
                if 'MockMC' in sample.shortName:
                    h_temp = files[channel].Get(f"h_SymmBkg{channel_charge}{pt_bin}_postFit")
                    if not h_temp:
                        btag = re.findall("([012]tag)", channel.name)[0]
                        h_temp = files[channel].Get(f"h_SymmBkg{channel_charge}_{btag}{pt_bin}_postFit")
                else:
                    if "SS" in channel.name:
                        h_temp = files[channel].Get(f"h_{sample.shortName}{truth_pt_bin}{channel_charge}_postFit")
                        # h_temp = files[channel].Get(f"h_{sample.shortName}{channel_charge}_SS_postFit")
                        h_temp_couter = files[channel].Get(f"h_{sample.shortName}{truth_pt_bin}{channel_charge}_CounterTerm_SS_postFit")
                        if h_temp_couter:
                            h_temp.Add(h_temp_couter)
                    else:
                        h_temp = files[channel].Get(f"h_{sample.shortName}{truth_pt_bin}{channel_charge}_postFit")
                        h_temp_couter = files[channel].Get(f"h_{sample.shortName}{truth_pt_bin}{channel_charge}_CounterTerm_postFit")
                        if h_temp_couter:
                            print("subtractign the counter term: ", h_temp_couter)
                            h_temp.Add(h_temp_couter)
                if h_temp:
                    if h_sum is None:
                        h_sum = h_temp.Clone(f"{h_temp.GetName()}_{chan.name}")
                    else:
                        h_sum.Add(h_temp)

            # negative channels
            for channel in plot['-']:

                # channel charge
                channel_charge = ""
                if options.plus_minus_ratio and ('MockMC' in sample.shortName or 'Wjets_emu_Matched' in sample.shortName):
                    if "minus" in channel.name:
                        channel_charge = "_minus"
                    elif "plus" in channel.name:
                        channel_charge = "_plus"

                # pt bin
                pt_bin = ""
                if 'MockMC' in sample.shortName or 'Wjets_emu_Matched' in sample.shortName:
                    pt_bin = "_" + re.findall("(pt_bin[0-9])", channel.name)[0]

                # truth pt bin
                truth_pt_bin = ""
                if pt_bin and 'Wjets_emu_Matched' in sample.shortName and "truth_pt" not in sample.shortName:
                    truth_pt_bin = f"_truth{pt_bin}"

                # mockMC
                if 'MockMC' in sample.shortName:
                    h_temp = files[channel].Get(f"h_SymmBkg{channel_charge}{pt_bin}_postFit")
                    if not h_temp:
                        btag = re.findall("([012]tag)", channel.name)[0]
                        h_temp = files[channel].Get(f"h_SymmBkg{channel_charge}_{btag}{pt_bin}_postFit")
                else:
                    h_temp = files[channel].Get(f"h_{sample.shortName}{truth_pt_bin}{channel_charge}_postFit")
                    # h_temp = files[channel].Get(f"h_{sample.shortName}{channel_charge}_SS_postFit")
                    h_temp_couter = files[channel].Get(f"h_{sample.shortName}{truth_pt_bin}{channel_charge}_CounterTerm_SS_postFit")
                    if h_temp_couter:
                        print("subtractign the counter term: ", h_temp_couter)
                        h_temp.Add(h_temp_couter)
                if h_temp:
                    if h_sum is None:
                        h_sum = h_temp.Clone(f"{h_temp.GetName()}_{chan.name}")
                        h_sum.Scale(-1.)
                    else:
                        h_sum.Add(h_temp, -1)
            if h_sum and abs(h_sum.GetSum()) > 1e-2:
                mc_map[sample] = h_sum

        # get data
        h_data = None
        for channel in plot['+']:
            h_temp = files[channel].Get("h_Data")
            if not h_data:
                h_data = h_temp.Clone(f"{h_temp.GetName()}_{chan.name}")
            else:
                h_data.Add(h_temp)
            if subtract_manually and len(plot['-']):
                for channel_SS in channels:
                    if channel_SS.name == channel.name.replace("OS_", "SS_"):
                        h_temp_SS = files[channel_SS].Get("h_Data")
                        h_data.Add(h_temp_SS, -1)
        if not subtract_manually:
            for channel in plot['-']:
                h_temp = files[channel].Get("h_Data")
                h_data.Add(h_temp, -1)

        # get mc tot
        h_mc_tot = None
        for channel in plot['+']:
            h_temp = files[channel].Get("h_tot_postFit")
            if not h_mc_tot:
                h_mc_tot = h_temp.Clone(f"{h_temp.GetName()}_{chan.name}")
            else:
                h_mc_tot.Add(h_temp)
            if subtract_manually and len(plot['-']):
                for channel_SS in channels:
                    if channel_SS.name == channel.name.replace("OS_", "SS_"):
                        h_temp_SS = files[channel_SS].Get("h_tot_postFit")
                        h_mc_tot.Add(h_temp_SS, -1)
        if not subtract_manually:
            for channel in plot['-']:
                h_temp = files[channel].Get("h_tot_postFit")
                h_mc_tot.Add(h_temp, -1)

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
            print(f"{par} {h_mc_tot.GetSumOfWeights()}")
            for channel in plot['+']:
                h_temp_up = get_err_hist(files[channel], par, "Up", "h_tot_postFit")
                h_temp_dn = get_err_hist(files[channel], par, "Down", "h_tot_postFit")
                if h_temp_up:
                    # print (f"plus up {h_temp_up.GetSumOfWeights()}")
                    if not h_mc_tot_Up:
                        h_mc_tot_Up = h_temp_up.Clone(f"{h_temp_up.GetName()}_{chan.name}_err_up")
                    else:
                        h_mc_tot_Up.Add(h_temp_up)
                if h_temp_dn:
                    # print (f"plus dn {h_temp_dn.GetSumOfWeights()}")
                    if not h_mc_tot_Dn:
                        h_mc_tot_Dn = h_temp_dn.Clone(f"{h_temp_dn.GetName()}_{chan.name}_err_dn")
                    else:
                        h_mc_tot_Dn.Add(h_temp_dn)
            if subtract_manually:
                for channel in plot['-']:
                    h_temp_up = get_err_hist(files[channel], par, "Up", "h_tot_postFit")
                    h_temp_dn = get_err_hist(files[channel], par, "Down", "h_tot_postFit")
                    if h_temp_up:
                        # print (f"minus up {h_temp_up.GetSumOfWeights()}")
                        if not h_mc_tot_Up:
                            h_mc_tot_Up = h_temp_up.Clone(f"{h_temp_up.GetName()}_{chan.name}_err_up")
                            h_mc_tot_Up.Scale(-1.)
                        else:
                            h_mc_tot_Up.Add(h_temp_up, -1)
                    if h_temp_dn:
                        # print (f"minus dn {h_temp_dn.GetSumOfWeights()}")
                        if not h_mc_tot_Dn:
                            h_mc_tot_Dn = h_temp_dn.Clone(f"{h_temp_dn.GetName()}_{chan.name}_err_dn")
                            h_mc_tot_Dn.Scale(-1)
                        else:
                            h_mc_tot_Dn.Add(h_temp_dn, -1)

            # subtract nominal
            h_mc_tot_Up.Add(h_mc_tot, -1)
            h_mc_tot_Dn.Add(h_mc_tot, -1)
            h_mc_tot_err_histograms_Up += [h_mc_tot_Up]
            h_mc_tot_err_histograms_Dn += [h_mc_tot_Dn]

        # calculate error bands
        for x in range(1, h_mc_tot.GetNbinsX() + 1):
            g_mc_tot_err.GetEYlow()[x - 1] = 0.
            g_mc_tot_err.GetEYhigh()[x - 1] = 0.
            # off diagonal
            for i in range(n_pars):
                for j in range(i):
                    if h_mc_tot_err_histograms_Up[i] and h_mc_tot_err_histograms_Up[j] and h_mc_tot_err_histograms_Dn[i] and h_mc_tot_err_histograms_Dn[j]:
                        corr = corr_correlation_rows[i][j]
                        # err_i = (h_mc_tot_err_histograms_Up[i].GetBinContent(x) - h_mc_tot_err_histograms_Dn[i].GetBinContent(x)) / 2.
                        # err_j = (h_mc_tot_err_histograms_Up[j].GetBinContent(x) - h_mc_tot_err_histograms_Dn[j].GetBinContent(x)) / 2.
                        # err = err_i * err_j * corr * 2
                        # g_mc_tot_err.GetEYlow()[x - 1] += err
                        # g_mc_tot_err.GetEYhigh()[x - 1] += err
                        err_up_i = h_mc_tot_err_histograms_Up[i].GetBinContent(x)
                        err_dn_i = h_mc_tot_err_histograms_Dn[i].GetBinContent(x)
                        err_up_j = h_mc_tot_err_histograms_Up[j].GetBinContent(x)
                        err_dn_j = h_mc_tot_err_histograms_Dn[j].GetBinContent(x)
                        err_up = err_up_i * err_up_j * corr * 2
                        err_dn = err_dn_i * err_dn_j * corr * 2
                        g_mc_tot_err.GetEYlow()[x - 1] += err_up
                        g_mc_tot_err.GetEYhigh()[x - 1] += err_dn
            # diagonal
            for i in range(n_pars):
                if h_mc_tot_err_histograms_Up[i] and h_mc_tot_err_histograms_Dn[i]:
                    # err_i = (h_mc_tot_err_histograms_Up[i].GetBinContent(x) - h_mc_tot_err_histograms_Dn[i].GetBinContent(x)) / 2.
                    # err = err_i * err_i
                    # g_mc_tot_err.GetEYlow()[x - 1] += err
                    # g_mc_tot_err.GetEYhigh()[x - 1] += err
                    err_up_i = h_mc_tot_err_histograms_Up[i].GetBinContent(x)
                    err_dn_i = h_mc_tot_err_histograms_Dn[i].GetBinContent(x)
                    err_up = err_up_i * err_up_i
                    err_dn = err_dn_i * err_dn_i
                    g_mc_tot_err.GetEYlow()[x - 1] += err_up
                    g_mc_tot_err.GetEYhigh()[x - 1] += err_dn

            # final
            g_mc_tot_err.GetEYhigh()[x - 1] = math.sqrt(g_mc_tot_err.GetEYhigh()[x - 1])
            g_mc_tot_err.GetEYlow()[x - 1] = math.sqrt(g_mc_tot_err.GetEYlow()[x - 1])
            g_mc_tot_err_only.GetEYhigh()[x - 1] = g_mc_tot_err.GetEYhigh()[x - 1] / h_mc_tot.GetBinContent(x)
            g_mc_tot_err_only.GetEYlow()[x - 1] = g_mc_tot_err.GetEYlow()[x - 1] / h_mc_tot.GetBinContent(x)

        # ratio
        h_ratio = utils.make_ratio(h_data, h_mc_tot)

        # data graph
        gr_data = utils.get_gr_from_hist(h_data)
        gr_ratio = utils.get_gr_from_hist(h_ratio)

        # top pad
        canv.pad1.cd()

        # normalize bins to unity
        if var.per_unit:
            utils.normalize_to_unit(hs, hists=[h_mc_tot], grs=[g_mc_tot_err, gr_data])

        hs.Draw("same hist")
        h_mc_tot.Draw("same hist")
        g_mc_tot_err.Draw("e2")
        gr_data.Draw("pe0")

        # make legend
        canv.make_legend(h_data, h_mc_tot, mc_map, samples, print_yields=True, show_error=False)

        # set maximum after creating legend
        canv.set_maximum((h_data, h_mc_tot), var, mc_min=utils.get_mc_min(mc_map, samples))

        # find minimum
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
        canv.set_ratio_range(0.51, 1.49, override=True)
        ROOT.gPad.RedrawAxis()

        # Print out
        canv.print(f"{conf.out_name}/{chan.name}_{var.name}.pdf")
        canv.print(f"{conf.out_name}/{chan.name}_{var.name}.png")

        logging.info(f"finished processing channel {channel.name}")


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-a', '--analysis-config',
                      action="store", dest="analysis_config",
                      help="analysis config file")
    parser.add_option('-v', '--var',
                      action="store", dest="var",
                      help="fitted variable",
                      default="Dmeson_m")
    parser.add_option('--suffix',
                      action="store", dest="suffix",
                      help="suffix for the output name")
    parser.add_option('--stage-out',
                      action="store_true", dest="stage_out",
                      help="copy plots to the www folder")
    parser.add_option('--plus-minus-ratio',
                      action="store_true", dest="plus_minus_ratio",
                      help="the fit was a combined plus minus fit")
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
    if not os.path.isdir(out_name):
        os.makedirs(out_name)

    # do the plotting
    main(options, conf)

    # stage-out to the www folder
    if options.stage_out:
        www.stage_out_plots(conf.out_name, conf.get_variables(), x=300, y=300)
