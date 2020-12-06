#!/usr/bin/env python
from charmplot.control import globalConfig
from charmplot.control import tools
from charmplot.control.channel import Channel
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
        h_mc_default = f.Get(default)
        if h_mc_default:
            h_err = h_mc_default.Clone(f"{h_mc_default.GetName()}_err_{variation}")
    return h_err


def main(options, conf):
    trex_histogram_folder = os.path.join(options.trex_input, "Histograms")

    # channels
    channels = []

    # parse input
    for file in os.listdir(trex_histogram_folder):
        if file.endswith("postFit.root"):
            channel_name = file.replace("_postFit.root", "")
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
    for channel_OS in channels:
        if "OS_" not in channel_OS.name:
            continue
        print(channel_OS.name)
        for channel_SS in channels:
            if channel_SS.name == channel_OS.name.replace("OS_", "SS_"):
                OS_minus_SS_plots += [{'+': [channel_OS], '-': [channel_SS]}]
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

    # plots = individual_plots + OS_minus_SS_plots + [OS_minus_SS_total]
    plots = OS_minus_SS_plots

    for plot in plots:

        # create channel
        channel_temp = plot['+'][0]
        channels_all = plot['+'] + plot['-']
        channel_name = "_".join([channel.name for channel in channels_all])
        chan = Channel(channel_name, [], channel_temp.lumi, [], [])

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
                    if options.samples and sample.shortName not in options.samples.split(","):
                        continue
                    sample_names.append(sample.shortName)
                    samples.append(sample)

        # get mc samples
        mc_map = {}
        sample_names = {}
        for sample in samples:
            sample_names[sample] = {}
            h_sum = None
            for channel in plot['+']:
                if 'MockMC' in sample.shortName:
                    btag = re.findall("([012]tag)", channel.name)[0]
                    name = f"h_{sample.shortName}_{btag}_postFit"
                    h_temp = files[channel].Get(name)
                    sample_names[sample][channel] = name
                else:
                    if "SS" in channel.name:
                        name = f"h_{sample.shortName}_SS_postFit"
                        h_temp = files[channel].Get(name)
                    else:
                        name = f"h_{sample.shortName}_postFit"
                        h_temp = files[channel].Get(name)
                    sample_names[sample][channel] = name
                if h_temp:
                    if h_sum is None:
                        h_sum = h_temp.Clone(f"{h_temp.GetName()}_{chan.name}")
                    else:
                        h_sum.Add(h_temp)
            for channel in plot['-']:
                if 'MockMC' in sample.shortName:
                    btag = re.findall("([012]tag)", channel.name)[0]
                    name = f"h_{sample.shortName}_{btag}_postFit"
                    h_temp = files[channel].Get(name)
                else:
                    name = f"h_{sample.shortName}_SS_postFit"
                    h_temp = files[channel].Get(name)
                sample_names[sample][channel] = name
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
        for channel in plot['-']:
            h_temp = files[channel].Get("h_tot_postFit")
            h_mc_tot.Add(h_temp, -1)

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
            h_mc_tot_Up.Add(h_mc_tot, -1)
            h_mc_tot_Dn.Add(h_mc_tot, -1)
            h_mc_tot_err_histograms_Up += [h_mc_tot_Up]
            h_mc_tot_err_histograms_Dn += [h_mc_tot_Dn]

        # calculate error bands
        h_mc_tot_Err = 0
        # off diagonal
        for i in range(n_pars):
            for j in range(i):
                if h_mc_tot_err_histograms_Up[i] and h_mc_tot_err_histograms_Up[j] and h_mc_tot_err_histograms_Dn[i] and h_mc_tot_err_histograms_Dn[j]:
                    corr = corr_correlation_rows[i][j]
                    err_i = (h_mc_tot_err_histograms_Up[i].GetSum() - h_mc_tot_err_histograms_Dn[i].GetSum()) / 2.
                    err_j = (h_mc_tot_err_histograms_Up[j].GetSum() - h_mc_tot_err_histograms_Dn[j].GetSum()) / 2.
                    err = err_i * err_j * corr * 2
                    h_mc_tot_Err += err
        # diagonal
        for i in range(n_pars):
            if h_mc_tot_err_histograms_Up[i] and h_mc_tot_err_histograms_Dn[i]:
                err_i = (h_mc_tot_err_histograms_Up[i].GetSum() - h_mc_tot_err_histograms_Dn[i].GetSum()) / 2.
                err = err_i * err_i
                h_mc_tot_Err += err

        # final
        h_mc_tot_Err = math.sqrt(h_mc_tot_Err)

        # uncertainties for each sample
        sample_errors = {}
        for sample in samples:
            # systematic uncertainties
            h_sample_err_histograms_Up = []
            h_sample_err_histograms_Dn = []
            for par in corr_parameters:
                h_sample_Up = None
                h_sample_Dn = None
                for channel in plot['+']:
                    h_temp_up = get_err_hist(files[channel], par, "Up", sample_names[sample][channel])
                    h_temp_dn = get_err_hist(files[channel], par, "Down", sample_names[sample][channel])
                    if h_temp_up:
                        if not h_sample_Up:
                            h_sample_Up = h_temp_up.Clone(f"{h_temp_up.GetName()}_{chan.name}_err_up")
                        else:
                            h_sample_Up.Add(h_temp_up)
                    if h_temp_dn:
                        if not h_sample_Dn:
                            h_sample_Dn = h_temp_dn.Clone(f"{h_temp_dn.GetName()}_{chan.name}_err_dn")
                        else:
                            h_sample_Dn.Add(h_temp_dn)
                for channel in plot['-']:
                    h_temp_up = get_err_hist(files[channel], par, "Up", sample_names[sample][channel])
                    h_temp_dn = get_err_hist(files[channel], par, "Down", sample_names[sample][channel])
                    if h_temp_up:
                        if not h_sample_Up:
                            h_sample_Up = h_temp_up.Clone(f"{h_temp_up.GetName()}_{chan.name}_err_up")
                            h_sample_Up.Scale(-1.)
                        else:
                            h_sample_Up.Add(h_temp_up, -1)
                    if h_temp_dn:
                        if not h_sample_Dn:
                            h_sample_Dn = h_temp_dn.Clone(f"{h_temp_dn.GetName()}_{chan.name}_err_dn")
                        else:
                            h_sample_Dn.Add(h_temp_dn, -1)
                            h_sample_Dn.Scale(-1.)

                # subtract nominal
                if sample in mc_map:
                    h_sample_Up.Add(mc_map[sample], -1)
                    h_sample_Dn.Add(mc_map[sample], -1)
                h_sample_err_histograms_Up += [h_sample_Up]
                h_sample_err_histograms_Dn += [h_sample_Dn]

            # calculate error bands
            h_sample_Err = 0
            h_sample_Err_matrix = []
            # off diagonal
            for i in range(n_pars):
                for j in range(i):
                    if h_sample_err_histograms_Up[i] and h_sample_err_histograms_Up[j] and h_sample_err_histograms_Dn[i] and h_sample_err_histograms_Dn[j]:
                        corr = corr_correlation_rows[i][j]
                        err_i = (h_sample_err_histograms_Up[i].GetSum() - h_sample_err_histograms_Dn[i].GetSum()) / 2.
                        err_j = (h_sample_err_histograms_Up[j].GetSum() - h_sample_err_histograms_Dn[j].GetSum()) / 2.
                        err = err_i * err_j * corr * 2
                        h_sample_Err += err
                        h_sample_Err_matrix += [[err, (i, j)]]
            # diagonal
            for i in range(n_pars):
                if h_sample_err_histograms_Up[i] and h_sample_err_histograms_Dn[i]:
                    err_i = (h_sample_err_histograms_Up[i].GetSum() - h_sample_err_histograms_Dn[i].GetSum()) / 2.
                    err = err_i * err_i
                    h_sample_Err += err
                    h_sample_Err_matrix += [[err, (i, i)]]

            # final
            h_sample_Err = math.sqrt(h_sample_Err)
            h_sample_Err_matrix = sorted(h_sample_Err_matrix, key=lambda x: x[0], reverse=True)
            sample_errors[sample] = h_sample_Err

        # print out
        logging.info(f"=== Printing yields for {chan.name} ===")
        logging.info(f"Data: {h_data.GetSum()}")
        logging.info(f"MC. Tot.: {h_mc_tot.GetSum()} += {h_mc_tot_Err}")
        for sample in samples:
            if sample in mc_map and sample in sample_errors:
                logging.info(f"{sample.shortName}: {mc_map[sample].GetSum()} += {sample_errors[sample]}")


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-a', '--analysis-config',
                      action="store", dest="analysis_config",
                      help="analysis config file")
    parser.add_option('-s', '--samples',
                      action="store", dest="samples",
                      help="get yields only for the specified samples")
    parser.add_option('--trex-input',
                      action="store", dest="trex_input",
                      help="import post-fit trex plots")

    # parse input arguments
    options, args = parser.parse_args()

    # analysis configs
    config = options.analysis_config

    # output name
    out_name = config.replace(".yaml", "").replace(".yml", "")

    # config object
    conf = globalConfig.GlobalConfig(config, out_name)

    # do the plotting
    main(options, conf)
