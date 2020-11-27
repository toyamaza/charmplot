#!/usr/bin/env python
from charmplot.common import utils
from charmplot.common import www
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


def get_err_hist(f, par):
    h_err_up = None
    h_temp = f.Get(f"h_tot_{par}_Up_postFit")
    if h_temp:
        h_err_up = h_temp.Clone(f"{h_temp.GetName()}_err")
        h_err_up.Add(f.Get("h_tot_postFit"), -1)
    return h_err_up


def main(options, conf, total=False):
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
            channel = conf.get_channel(channel_name)
            if not channel:
                logging.error(f"Channel not found for string {channel_name}")
            else:
                logging.info(f"Found channel {channel_name}")
            channels += [channel]

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

    # loop over channels
    for channel in channels:

        # only OS
        if not channel.name.startswith("OS_"):
            continue

        if total and "0tag" not in channel.name:
            continue

        # labels
        if not total:
            labels = []
            for lab in channel.label:
                lab = lab.replace("OS", "OS-SS")
                labels += [lab]
            channel.label = labels + ["post-fit"]
        else:
            labels = []
            for lab in channel.label:
                lab = lab.replace(", 0tag", "")
                labels += [lab]
            channel.label = labels

        # samples
        samples = []

        # read root file
        f = ROOT.TFile(os.path.join(trex_histogram_folder, f"{channel.name}_postFit.root"))
        f_SS = ROOT.TFile(os.path.join(trex_histogram_folder, f"{channel.name.replace('OS_', 'SS_')}_postFit.root"))
        if total:
            f_1 = ROOT.TFile(os.path.join(trex_histogram_folder, f"{channel.name.replace('0tag', '1tag')}_postFit.root"))
            f_SS_1 = ROOT.TFile(os.path.join(trex_histogram_folder, f"{channel.name.replace('OS_', 'SS_').replace('0tag', '1tag')}_postFit.root"))
            f_2 = ROOT.TFile(os.path.join(trex_histogram_folder, f"{channel.name.replace('0tag', '2tag')}_postFit.root"))
            f_SS_2 = ROOT.TFile(os.path.join(trex_histogram_folder, f"{channel.name.replace('OS_', 'SS_').replace('0tag', '2tag')}_postFit.root"))

        # get mc samples
        mc_map = {}
        for sample_name in channel.samples:
            if "MockMC" in sample_name:
                continue
            sample = conf.get_sample(sample_name)
            samples += [sample]
            h_name = f"h_{sample.shortName}_postFit"
            h = f.Get(h_name)
            if not h:
                continue
            h_SS = f_SS.Get(h_name)
            if h_SS:
                h.Add(h_SS, -1)
            logging.info(f"Got histogram {h} for sample {sample_name}")
            if total:
                h.Add(f_1.Get(h_name))
                h.Add(f_2.Get(h_name))
                if f_SS_1.Get(h_name):
                    h.Add(f_SS_1.Get(h_name), -1)
                if f_SS_2.Get(h_name):
                    h.Add(f_SS_2.Get(h_name), -1)
            mc_map[sample] = h

        # get data
        h_data = f.Get("h_Data")
        h_data_SS = f_SS.Get("h_Data")
        h_data.Add(h_data_SS, -1)
        if total:
            h_data.Add(f_1.Get("h_Data"))
            h_data.Add(f_2.Get("h_Data"))
            h_data.Add(f_SS_1.Get("h_Data"), -1)
            h_data.Add(f_SS_2.Get("h_Data"), -1)
        logging.info(f"Got histogram {h_data} for Data")

        # canvas
        canv = utils.make_canvas(h_data, var, channel, x=800, y=800)

        # configure histograms
        canv.configure_histograms(mc_map, h_data)

        # stack and total mc
        hs = utils.make_stack(samples, mc_map)
        h_mc_tot_OS = f.Get("h_tot_postFit")
        h_mc_tot_SS = f_SS.Get("h_tot_postFit")
        h_mc_tot = h_mc_tot_OS.Clone(f"{channel.name}_{var.name}_mc_tot")
        h_mc_tot.Add(h_mc_tot_SS, -1)
        h_mc_tot.SetLineWidth(1)
        h_mc_tot.SetFillStyle(0)
        h_mc_tot.SetLineColor(ROOT.kWhite)
        if total:
            h_mc_tot.Add(f_1.Get("h_tot_postFit"))
            h_mc_tot.Add(f_2.Get("h_tot_postFit"))
            h_mc_tot.Add(f_SS_1.Get("h_tot_postFit"), -1)
            h_mc_tot.Add(f_SS_2.Get("h_tot_postFit"), -1)

        # stat error
        g_mc_tot_err, g_mc_tot_err_only = utils.make_stat_err(h_mc_tot)

        # systematic uncertainties
        h_mc_tot_err_histograms = []
        for par in corr_parameters:
            h_err = None
            h_sys_OS_up = get_err_hist(f, par)
            h_sys_SS_up = get_err_hist(f_SS, par)
            if h_sys_OS_up:
                h_err = h_sys_OS_up.Clone(f"{h_sys_OS_up.GetName()}_final")
            if h_sys_SS_up and h_err:
                h_err.Add(h_sys_SS_up, -1)
            if total:
                # 1-tag
                h_sys_OS_up_1 = get_err_hist(f_1, par)
                h_sys_SS_up_1 = get_err_hist(f_SS_1, par)
                if h_sys_OS_up_1:
                    if h_err:
                        h_err.Add(h_sys_OS_up_1)
                    else:
                        h_err = h_sys_OS_up_1.Clone(f"{h_sys_OS_up_1.GetName()}_final")
                if h_sys_SS_up_1 and h_err:
                    h_err.Add(h_sys_SS_up_1, -1)
                # 2-tag
                h_sys_OS_up_2 = get_err_hist(f_2, par)
                h_sys_SS_up_2 = get_err_hist(f_SS_2, par)
                if h_sys_OS_up_2:
                    if h_err:
                        h_err.Add(h_sys_OS_up_2)
                    else:
                        h_err = h_sys_OS_up_2.Clone(f"{h_sys_OS_up_2.GetName()}_final")
                if h_sys_SS_up_2 and h_err:
                    h_err.Add(h_sys_SS_up_2, -1)

            h_mc_tot_err_histograms += [h_err]

        for x in range(1, h_mc_tot.GetNbinsX() + 1):
            g_mc_tot_err.GetEYlow()[x] = 0.
            g_mc_tot_err.GetEYhigh()[x] = 0.
            # off diagonal
            for i in range(n_pars):
                for j in range(i):
                    if h_mc_tot_err_histograms[i] and h_mc_tot_err_histograms[j]:
                        corr = corr_correlation_rows[i][j]
                        err_i = h_mc_tot_err_histograms[i].GetBinContent(x)
                        err_j = h_mc_tot_err_histograms[j].GetBinContent(x)
                        err = err_i * err_j * corr * 2
                        g_mc_tot_err.GetEYlow()[x] += err
                        g_mc_tot_err.GetEYhigh()[x] += err
            # diagonal
            for i in range(n_pars):
                if h_mc_tot_err_histograms[i]:
                    err_i = h_mc_tot_err_histograms[i].GetBinContent(x)
                    err = err_i * err_i
                    g_mc_tot_err.GetEYlow()[x] += err
                    g_mc_tot_err.GetEYhigh()[x] += err

            # final
            g_mc_tot_err.GetEYhigh()[x] = math.sqrt(g_mc_tot_err.GetEYhigh()[x])
            g_mc_tot_err.GetEYlow()[x] = math.sqrt(g_mc_tot_err.GetEYlow()[x])
            g_mc_tot_err_only.GetEYhigh()[x] = g_mc_tot_err.GetEYhigh()[x] / h_mc_tot.GetBinContent(x)
            g_mc_tot_err_only.GetEYlow()[x] = g_mc_tot_err.GetEYlow()[x] / h_mc_tot.GetBinContent(x)

        # ratio
        h_ratio = utils.make_ratio(h_data, h_mc_tot)

        # top pad
        canv.pad1.cd()
        hs.Draw("same hist")
        h_mc_tot.Draw("same hist")
        g_mc_tot_err.Draw("e2")
        h_data.Draw("same pe")

        # make legend
        canv.make_legend(h_data, h_mc_tot, mc_map, samples, print_yields=True)

        # set maximum after creating legend
        canv.set_maximum((h_data, h_mc_tot), var, mc_min=utils.get_mc_min(mc_map, samples))

        # find minimum
        min_negative = {}
        for s in samples:
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
        canv.pad2.cd()
        g_mc_tot_err_only.Draw("le2")
        h_ratio.Draw("same pe")

        # Print out
        if not total:
            canv.print(f"{conf.out_name}/{channel.name}_{var.name}.pdf")
        else:
            canv.print(f"{conf.out_name}/{channel.name}_{var.name}_total.pdf")

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
    main(options, conf, total=True)

    # stage-out to the www folder
    if options.stage_out:
        www.stage_out_plots(conf.out_name, conf.get_variables(), x=300, y=300)
