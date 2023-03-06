#!/usr/bin/env python
from charmplot.control import globalConfig
from charmplot.control import tools
from charmplot.control.channel import Channel
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


def modify_sample_name(sample_name, channel):
    return sample_name


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


def get_yield_error(plot, chan, corr_parameters, corr_correlation_rows, result, files, h_ref, h_name="h_tot"):
    n_pars = len(corr_parameters)
    h_err_histograms_Up = []
    h_err_histograms_Dn = []
    for par in corr_parameters:
        h_Up = None
        h_Dn = None
        for channel in plot['+']:
            h_temp_up = get_err_hist(files[channel], par, "Up", f"{h_name}_postFit")
            h_temp_dn = get_err_hist(files[channel], par, "Down", f"{h_name}_postFit")
            if h_temp_up:
                if not h_Up:
                    h_Up = h_temp_up.Clone(f"{h_temp_up.GetName()}_{chan.name}_err_up")
                else:
                    h_Up.Add(h_temp_up)
            if h_temp_dn:
                if not h_Dn:
                    h_Dn = h_temp_dn.Clone(f"{h_temp_dn.GetName()}_{chan.name}_err_dn")
                else:
                    h_Dn.Add(h_temp_dn)
        for channel in plot['-']:
            h_temp_up = get_err_hist(files[channel], par, "Up", f"{h_name}_postFit")
            h_temp_dn = get_err_hist(files[channel], par, "Down", f"{h_name}_postFit")
            if h_temp_up:
                if not h_Up:
                    h_Up = h_temp_up.Clone(f"{h_temp_up.GetName()}_{chan.name}_err_up")
                    h_Up.Scale(-1.)
                else:
                    h_Up.Add(h_temp_up, -1)
            if h_temp_dn:
                if not h_Dn:
                    h_Dn = h_temp_dn.Clone(f"{h_temp_dn.GetName()}_{chan.name}_err_dn")
                    h_Dn.Scale(-1)
                else:
                    h_Dn.Add(h_temp_dn, -1)

        # subtract nominal
        h_Up.Add(h_ref, -1)
        h_Dn.Add(h_ref, -1)
        h_err_histograms_Up += [h_Up]
        h_err_histograms_Dn += [h_Dn]

    errors = []
    mc_tot_err_low = 0.
    mc_tot_err_high = 0.
    # off diagonal
    for i in range(n_pars):
        for j in range(i):
            # corr = corr_correlation_rows[i][j]
            par1 = corr_parameters[i]
            par2 = corr_parameters[j]
            if "SymmBkg" not in par1 and "mu_" not in par1 and not par1.startswith("stat_"):
                par1 = "alpha_" + par1
            elif "SymmBkg" in par1 or par1.startswith("stat_"):
                par1 = "gamma_" + par1
            if "SymmBkg" not in par2 and "mu_" not in par2 and not par2.startswith("stat_"):
                par2 = "alpha_" + par2
            elif "SymmBkg" in par2 or par2.startswith("stat_"):
                par2 = "gamma_" + par2
            corr = result.correlation(par1, par2)
            err_up_i = 0.
            err_dn_i = 0.
            err_up_j = 0.
            err_dn_j = 0.
            if h_err_histograms_Up[i]:
                err_up_i = h_err_histograms_Up[i].GetSumOfWeights()
            if h_err_histograms_Dn[i]:
                err_dn_i = h_err_histograms_Dn[i].GetSumOfWeights()
            if h_err_histograms_Up[j]:
                err_up_j = h_err_histograms_Up[j].GetSumOfWeights()
            if h_err_histograms_Dn[j]:
                err_dn_j = h_err_histograms_Dn[j].GetSumOfWeights()
            err_up_i = (err_up_i - err_dn_i) / 2.
            err_dn_i = err_up_i
            err_up_j = (err_up_j - err_dn_j) / 2.
            err_dn_j = err_up_j
            err_up = err_up_i * err_up_j * corr * 2.
            err_dn = err_dn_i * err_dn_j * corr * 2.
            mc_tot_err_low += err_dn
            mc_tot_err_high += err_up
            errors += [[err_up, (corr_parameters[i], corr_parameters[j], err_up_i, err_up_j, corr)]]
    # diagonal
    for i in range(n_pars):
        err_up_i = 0.
        err_dn_i = 0.
        if h_err_histograms_Up[i]:
            err_up_i = h_err_histograms_Up[i].GetSumOfWeights()
        if h_err_histograms_Dn[i]:
            err_dn_i = h_err_histograms_Dn[i].GetSumOfWeights()
        err_up_i = (err_up_i - err_dn_i) / 2.
        err_dn_i = err_up_i
        err_up = err_up_i * err_up_i
        err_dn = err_dn_i * err_dn_i
        mc_tot_err_low += err_dn
        mc_tot_err_high += err_up
        errors += [[err_up, (corr_parameters[i], corr_parameters[i], 0.0)]]

    # final
    err = c_double(0)
    _ = h_ref.IntegralAndError(0, h_ref.GetNbinsX() + 1, err)
    mc_tot_err_low = math.sqrt(mc_tot_err_low + err.value * err.value)
    mc_tot_err_high = math.sqrt(mc_tot_err_high + err.value * err.value)
    errors.sort(key=lambda x: abs(x[0]), reverse=True)
    return mc_tot_err_low, mc_tot_err_high, errors


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
                logging.warning(f"Channel not found for string {channel_name}")
            else:
                logging.info(f"Found channel {channel_name}")
                if options.skip_channel and (options.skip_channel in channel.name):
                    continue
                else:
                    channels += [channel]

    # sort channels
    individual_plots = []
    OS_minus_SS_plots = []
    OS_minus_SS_total = {'+': [], '-': []}
    OS_total = {'+': [], '-': []}
    SS_total = {'+': [], '-': []}
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
                    OS_total['+'] += [channel_OS]
                    SS_total['+'] += [channel_SS]
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
    res_file = ROOT.TFile(os.path.join(options.trex_input, "Fits", os.path.basename(options.trex_input) + ".root"), "READ")
    result = res_file.Get("nll_simPdf_newasimovData_with_constr")
    assert result
    corr_parameters = None
    corr_correlation_rows = None
    for x in corr_dict:
        if 'parameters' in x:
            corr_parameters = x['parameters']
        elif 'correlation_rows' in x:
            corr_correlation_rows = x['correlation_rows']

    # plots = individual_plots + OS_minus_SS_plots
    # plots = OS_minus_SS_plots + [OS_minus_SS_total]
    plots = OS_minus_SS_plots + [OS_minus_SS_total_minus, OS_minus_SS_total_plus]
    # plots = [OS_minus_SS_total, OS_total, SS_total]
    # plots = OS_minus_SS_plots + [OS_minus_SS_total, OS_total, SS_total]
    # plots = individual_plots + OS_minus_SS_plots + [OS_minus_SS_total, OS_total, SS_total]
    # plots = individual_plots + OS_minus_SS_plots + [OS_minus_SS_total, OS_minus_SS_total_minus, OS_minus_SS_total_plus]
    # plots = [OS_minus_SS_total]

    # sort regions
    plots.sort(key=lambda x: x['+'][0].name)
    for x in plots:
        print(x)

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
        has_OS = False
        has_SS = False
        if len(channels_all) > 2:
            for c in channels_all:
                if "minus" in c.name:
                    minus = True
                elif "plus" in c.name:
                    plus = True
                if "OS" in c.name:
                    has_OS = True
                elif "SS" in c.name:
                    has_SS = True
            if minus and plus:
                inclusive = True
            elif minus and not plus:
                inclusive_minus = True
            elif plus and not minus:
                inclusive_plus = True

        if inclusive:
            channel_name = "0tag_inclusive"
        elif inclusive_minus:
            channel_name = "0tag_inclusive_minus"
        elif inclusive_plus:
            channel_name = "0tag_inclusive_plus"

        if has_OS and not has_SS:
            channel_name = "OS_" + channel_name
        elif has_SS and not has_OS:
            channel_name = "SS_" + channel_name
        elif has_OS and has_SS:
            channel_name = "OS_SS_" + channel_name

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
                    sample_names.append(sample.shortName)
                    samples.append(sample)

        # get mc samples
        mc_map = {}
        for sample in samples:
            h_sum = None
            sample_name = sample.shortName

            # positive channels
            for channel in plot['+']:

                h_temp = files[channel].Get(f"h_{modify_sample_name(sample_name, channel)}_postFit")
                if h_temp:
                    if h_sum is None:
                        h_sum = h_temp.Clone(f"{h_temp.GetName()}_{chan.name}")
                    else:
                        h_sum.Add(h_temp)

            # negative channels
            for channel in plot['-']:
                h_temp = files[channel].Get(f"h_{modify_sample_name(sample_name, channel)}_postFit")
                if h_temp:
                    if h_sum is None:
                        h_sum = h_temp.Clone(f"{h_temp.GetName()}_{chan.name}")
                        h_sum.Scale(-1.)
                    else:
                        h_sum.Add(h_temp, -1)

            if h_sum and abs(h_sum.GetSum()) > 1e-6:
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

        # systematic uncertainties for MC tot
        mc_tot_err_low, mc_tot_err_high, _ = get_yield_error(plot, chan, corr_parameters, corr_correlation_rows, result, files, h_mc_tot, h_name="h_tot")

        # sys for signal
        print("Sys for signal 1")
        s_signal1 = [x for x in samples if x.shortName == f"Sherpa2211_WplusD_Matched_truth_{options.diff}_bin1"][0]
        if s_signal1 in mc_map:
            _, signal1_err_high, _ = get_yield_error(
                plot, chan, corr_parameters, corr_correlation_rows, result, files, mc_map[s_signal1], h_name=f"h_Sherpa2211_WplusD_Matched_truth_{options.diff}_bin1")
        else:
            signal1_err_high = 0

        # sys for signal
        print("Sys for signal 2")
        s_signal2 = [x for x in samples if x.shortName == f"Sherpa2211_WplusD_Matched_truth_{options.diff}_bin2"][0]
        if s_signal2 in mc_map:
            _, signal2_err_high, _ = get_yield_error(
                plot, chan, corr_parameters, corr_correlation_rows, result, files, mc_map[s_signal2], h_name=f"h_Sherpa2211_WplusD_Matched_truth_{options.diff}_bin2")
        else:
            signal2_err_high = 0

        # sys for signal
        print("Sys for signal 3")
        s_signal3 = [x for x in samples if x.shortName == f"Sherpa2211_WplusD_Matched_truth_{options.diff}_bin3"][0]
        if s_signal3 in mc_map:
            _, signal3_err_high, _ = get_yield_error(
                plot, chan, corr_parameters, corr_correlation_rows, result, files, mc_map[s_signal3], h_name=f"h_Sherpa2211_WplusD_Matched_truth_{options.diff}_bin3")
        else:
            signal3_err_high = 0

        # sys for signal
        print("Sys for signal 4")
        s_signal4 = [x for x in samples if x.shortName == f"Sherpa2211_WplusD_Matched_truth_{options.diff}_bin4"][0]
        if s_signal4 in mc_map:
            _, signal4_err_high, _ = get_yield_error(
                plot, chan, corr_parameters, corr_correlation_rows, result, files, mc_map[s_signal4], h_name=f"h_Sherpa2211_WplusD_Matched_truth_{options.diff}_bin4")
        else:
            signal4_err_high = 0

        # sys for signal
        print("Sys for signal 5")
        s_signal5 = [x for x in samples if x.shortName == f"Sherpa2211_WplusD_Matched_truth_{options.diff}_bin5"][0]
        if s_signal5 in mc_map:
            _, signal5_err_high, _ = get_yield_error(
                plot, chan, corr_parameters, corr_correlation_rows, result, files, mc_map[s_signal5], h_name=f"h_Sherpa2211_WplusD_Matched_truth_{options.diff}_bin5")
        else:
            signal5_err_high = 0

        # sys W+c(match)
        print("Sys for W+c(match)")
        if options.decay_mode == "Dplus":
            s_wcmatch = [x for x in samples if x.shortName == "MG_Wjets_Charm"][0]
            _, wcmatch_err_high, _ = get_yield_error(
                plot, chan, corr_parameters, corr_correlation_rows, result, files, mc_map[s_wcmatch], h_name="h_MG_Wjets_Charm")
        else:
            s_wcmatch = [x for x in samples if x.shortName == "Sherpa2211_Wjets_Charm"][0]
            _, wcmatch_err_high, _ = get_yield_error(
                plot, chan, corr_parameters, corr_correlation_rows, result, files, mc_map[s_wcmatch], h_name="h_Sherpa2211_Wjets_Charm")

        # sys W+c(mis-match)
        print("Sys for W+c(mis-match)")
        s_wcmismatch = [x for x in samples if x.shortName == "Sherpa2211_Wjets_MisMatched"][0]
        _, wcmismatch_err_high, _ = get_yield_error(
            plot, chan, corr_parameters, corr_correlation_rows, result, files, mc_map[s_wcmismatch], h_name="h_Sherpa2211_Wjets_MisMatched")

        # sys W+jets
        print("Sys for W+jets")
        if options.decay_mode == "Dplus":
            s_wjets = [x for x in samples if x.shortName == "Sherpa2211_Wjets_Rest"][0]
            _, wjets_err_high, _ = get_yield_error(plot, chan, corr_parameters, corr_correlation_rows, result,
                                                   files, mc_map[s_wjets], h_name="h_Sherpa2211_Wjets_Rest")
        else:
            s_wjets = [x for x in samples if x.shortName == "MG_Wjets_Rest"][0]
            _, wjets_err_high, _ = get_yield_error(plot, chan, corr_parameters, corr_correlation_rows, result,
                                                   files, mc_map[s_wjets], h_name="h_MG_Wjets_Rest")

        # sys Other
        print("Sys for Other")
        s_other = [x for x in samples if x.shortName == "DibosonZjets"][0]
        _, other_err_high, _ = get_yield_error(plot, chan, corr_parameters, corr_correlation_rows, result, files, mc_map[s_other], h_name="h_DibosonZjets")

        # sys Top
        print("Sys for Top")
        s_top = [x for x in samples if x.shortName == "Top"][0]
        _, top_err_high, _ = get_yield_error(plot, chan, corr_parameters, corr_correlation_rows, result, files, mc_map[s_top], h_name="h_Top")

        # sys for MultiJet
        print("Sys for MJ")
        s_mj = [x for x in samples if x.shortName == "Multijet_MatrixMethod"][0]
        _, multijet_err_high, _ = get_yield_error(
            plot, chan, corr_parameters, corr_correlation_rows, result, files, mc_map[s_mj], h_name="h_Multijet_MatrixMethod")

        # print
        print(f"--- {channel_name} ---")
        err = c_double(0)
        data_integral = h_data.IntegralAndError(0, h_data.GetNbinsX() + 1, err)
        if s_signal1 in mc_map:
            print(f"{'W+D(bin 1)':20s}: {mc_map[s_signal1].GetSumOfWeights():8.2f} & {signal1_err_high:8.2f}")
        if s_signal2 in mc_map:
            print(f"{'W+D(bin 2)':20s}: {mc_map[s_signal2].GetSumOfWeights():8.2f} & {signal2_err_high:8.2f}")
        if s_signal3 in mc_map:
            print(f"{'W+D(bin 3)':20s}: {mc_map[s_signal3].GetSumOfWeights():8.2f} & {signal3_err_high:8.2f}")
        if s_signal4 in mc_map:
            print(f"{'W+D(bin 4)':20s}: {mc_map[s_signal4].GetSumOfWeights():8.2f} & {signal4_err_high:8.2f}")
        if s_signal5 in mc_map:
            print(f"{'W+D(bin 5)':20s}: {mc_map[s_signal5].GetSumOfWeights():8.2f} & {signal5_err_high:8.2f}")
        print(f"{'W+c(match)':20s}: {mc_map[s_wcmatch].GetSumOfWeights():8.2f} & {wcmatch_err_high:8.2f}")
        print(f"{'W+c(mis-match)':20s}: {mc_map[s_wcmismatch].GetSumOfWeights():8.2f} & {wcmismatch_err_high:8.2f}")
        print(f"{'W+jets':20s}: {mc_map[s_wjets].GetSumOfWeights():8.2f} & {wjets_err_high:8.2f}")
        print(f"{'Top':20s}: {mc_map[s_top].GetSumOfWeights():8.2f} & {top_err_high:8.2f}")
        print(f"{'Other':20s}: {mc_map[s_other].GetSumOfWeights():8.2f} & {other_err_high:8.2f}")
        print(f"{'Multijet':20s}: {mc_map[s_mj].GetSumOfWeights():8.2f} & {multijet_err_high:8.2f}")
        print(f"{'SM tot.':20s}: {h_mc_tot.GetSumOfWeights():8.2f} & {mc_tot_err_high:8.2f}")
        print(f"{'Data':20s}: {data_integral:8.2f} & {err.value:8.2f}")

        for _, file in files.items():
            file.Close()
        logging.info(f"finished processing channel {channel_name}")


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
    parser.add_option('-d', '--diff',
                      action="store", dest="diff",
                      help="differential variable",
                      default="pt")
    parser.add_option('-k', '--skip-channel',
                      action="store", dest="skip_channel",
                      default="",
                      help="skip channels that include this in the name")
    parser.add_option('--trex-input',
                      action="store", dest="trex_input",
                      help="import post-fit trex plots")
    parser.add_option('--decay-mode',
                      action="store", dest="decay_mode",
                      default="Dplus")

    # parse input arguments
    options, args = parser.parse_args()

    # analysis configs
    config = options.analysis_config

    # output name
    out_name = config.replace(".yaml", "").replace(".yml", "")

    # config object
    conf = globalConfig.GlobalConfig(config, out_name)

    # make output folder if not exist
    if not os.path.isdir(os.path.join("post_fit", out_name)):
        os.makedirs(os.path.join("post_fit", out_name))

    # do the plotting
    main(options, conf)
