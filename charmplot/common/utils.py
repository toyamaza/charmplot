from charmplot.common import canvas
from charmplot.common import likelihoodFit
from charmplot.common import massFit
from charmplot.control import channel
from charmplot.control import globalConfig
from charmplot.control import inputDataReader
from charmplot.control import sample
from charmplot.control import variable
from copy import deepcopy
from ctypes import c_double
from typing import Dict, List, Union
import array
import lhapdf
import logging
import numpy as np
import os
import ROOT
import subprocess
import sys

MC_Map = Dict[sample.Sample, ROOT.TH1]

logger = logging.getLogger(__name__)


def hadd_wrapper(args):
    p = subprocess.Popen(["hadd"] + ["-f"] + args,
                         stderr=subprocess.STDOUT)
    p.communicate()
    return


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def get_gr_from_hist(h: ROOT.TH1) -> ROOT.TGraphErrors:
    gr = ROOT.TGraphErrors()
    for i in range(1, h.GetNbinsX() + 1):
        gr.SetPoint(gr.GetN(), h.GetBinCenter(i), h.GetBinContent(i))
        gr.SetPointError(gr.GetN() - 1, 0, h.GetBinError(i))
    gr.SetMarkerSize(1.0)
    return gr


def get_mc_min(mc_map: MC_Map, samples: list):
    for s in reversed(samples):
        if s not in mc_map.keys():
            continue
        return mc_map[s]


def read_sys_histograms(conf, reader, c, var, samples, fit, systematics, mc_map):
    mc_map_sys = {}
    for group in systematics:
        variations = systematics[group].get('variations')
        affecting = systematics[group].get('affecting')
        sys_type = systematics[group].get('type')
        if sys_type not in ['alt_sample', 'overall']:
            mc_map_sys[group] = {syst: read_samples(conf, reader, c, var, samples, fit,
                                                    force_positive=c.force_positive, sys=syst,
                                                    affecting=affecting, fallback=mc_map) for syst in variations}
        elif sys_type == 'alt_sample':
            # load files for alternative samples
            for syst in variations:
                sample = conf.construct_sample(syst)
                if not sample:
                    logging.error(f"sample not found for variation {syst}!")
                conf.add_sample(sample)
                reader.read_input_file(sample)
            # start with nominal for all samples
            mc_map_sys[group] = {syst: read_samples(conf, reader, c, var, samples, fit,
                                                    force_positive=c.force_positive, sys=syst,
                                                    affecting=[''], fallback=mc_map) for syst in variations}
            # replace affected with alternative samples
            for syst in variations:
                for s in affecting:
                    sample = next((x for x in samples if x.shortName == s), None)
                    if not sample:
                        continue
                    sample_sys = conf.get_sample(syst)
                    sample_sys.channel = sample.channel
                    h_sys = reader.get_histogram(sample_sys, c, var, c.force_positive)
                    if h_sys:
                        mc_map_sys[group][syst][sample] = h_sys
        elif sys_type == 'overall':
            # start with nominal for all samples
            mc_map_sys[group] = {syst: read_samples(conf, reader, c, var, samples, fit,
                                                    force_positive=c.force_positive, sys=syst,
                                                    affecting=[''], fallback=mc_map) for syst in variations}
            for syst in variations:
                for s in affecting:
                    sample = next((x for x in samples if x.shortName == s), None)
                    if not sample:
                        continue
                    if not sample in mc_map:
                        continue
                    h_nominal = mc_map[sample]
                    if not h_nominal:
                        continue
                    h_sys = h_nominal.Clone(f"{h_nominal.GetName()}_{syst}")
                    size = systematics[group].get('size')
                    if "1up" in syst:
                        h_sys.Scale(1 + size)
                    elif "1dn" in syst:
                        h_sys.Scale(1 - size)
                    mc_map_sys[group][syst][sample] = h_sys
    return mc_map_sys


def trex_subtraction(channelOS: channel.Channel, channelSS: channel.Channel,
                     var: variable.Variable, mc_map: MC_Map, trex_post_fit_histograms: Dict):
    logging.info(f"Making post-fit histograms for {channelOS.name} - {channelSS.name}")
    mc_map_OS = deepcopy(mc_map)
    mc_map_SS = deepcopy(mc_map)
    h_data_trex_OS, trex_mc_tot_OS, trex_mc_stat_err_OS, trex_mc_stat_err_only_OS = read_trex_input(channelOS, var, mc_map_OS, trex_post_fit_histograms)
    h_data_trex_SS, trex_mc_tot_SS, trex_mc_stat_err_SS, trex_mc_stat_err_only_SS = read_trex_input(channelSS, var, mc_map_SS, trex_post_fit_histograms)
    h_data_trex = h_data_trex_OS.Clone(h_data_trex_OS.GetName().replace("OS", "OS-SS"))
    h_data_trex.Add(h_data_trex_SS, -1)
    trex_mc_tot = trex_mc_tot_OS.Clone(trex_mc_tot_OS.GetName().replace("OS", "OS-SS"))
    trex_mc_tot.Add(trex_mc_tot_SS, -1)
    for s in mc_map_OS:
        h_temp = mc_map_OS[s].Clone(mc_map_OS[s].GetName().replace("OS", "OS-SS"))
        h_temp_SS = None
        for s2 in mc_map_SS:
            if s2.shortName == s.shortName:
                h_temp_SS = mc_map_SS[s2]
                break
        h_temp.Add(h_temp_SS, -1)
        logger.info(f"Subtracting {h_temp_SS.GetName()} {h_temp_SS.GetSum()} from {h_temp.GetSum()}")
        for s2 in mc_map:
            if s2.shortName == s.shortName:
                mc_map[s2] = h_temp
                break
    trex_mc_stat_err = trex_mc_stat_err_OS.Clone(trex_mc_stat_err_OS.GetName().replace("OS", "OS-SS"))
    trex_mc_stat_err_only = trex_mc_stat_err_only_OS.Clone(trex_mc_stat_err_only_OS.GetName().replace("OS", "OS-SS"))
    for i in range(trex_mc_stat_err.GetN()):
        yOS = trex_mc_stat_err_OS.GetY()[i]
        ySS = trex_mc_stat_err_SS.GetY()[i]
        yerr_lowOS = trex_mc_stat_err_OS.GetEYlow()[i]
        yerr_highOS = trex_mc_stat_err_OS.GetEYhigh()[i]
        yerr_lowSS = trex_mc_stat_err_SS.GetEYlow()[i]
        yerr_highSS = trex_mc_stat_err_SS.GetEYhigh()[i]
        trex_mc_stat_err.GetY()[i] = yOS - ySS
        trex_mc_stat_err.GetEYlow()[i] = yerr_lowOS + yerr_lowSS
        trex_mc_stat_err.GetEYhigh()[i] = yerr_highOS + yerr_highSS
        # trex_mc_stat_err.GetEYlow()[i] = (yerr_lowOS**2 + yerr_lowSS**2)**(0.5)
        # trex_mc_stat_err.GetEYhigh()[i] = (yerr_highOS**2 + yerr_highSS**2)**(0.5)
        if yOS - ySS > 0:
            trex_mc_stat_err_only.GetEYlow()[i] = trex_mc_stat_err.GetEYlow()[i] / (yOS - ySS)
            trex_mc_stat_err_only.GetEYhigh()[i] = trex_mc_stat_err.GetEYhigh()[i] / (yOS - ySS)
        else:
            trex_mc_stat_err_only.GetEYlow()[i] = 0
            trex_mc_stat_err_only.GetEYhigh()[i] = 0
    return h_data_trex, trex_mc_tot, trex_mc_stat_err, trex_mc_stat_err_only


def read_trex_input(channel: channel.Channel, var: variable.Variable, mc_map: MC_Map, trex_post_fit_histograms: Dict):
    file_temp = ROOT.TFile(trex_post_fit_histograms[channel.name], "READ")
    h_data_trex = deepcopy(file_temp.Get("h_Data"))
    var_trex = h_data_trex.GetTitle().split("__")[1]
    logging.info(f"TRExVar {var_trex}")
    if var.name == var_trex:
        logging.info(f"Replacing Data histogram with trex post-fit histogram {h_data_trex}")
        trex_mc_tot = deepcopy(file_temp.Get("h_tot_postFit"))
        trex_mc_tot.SetLineWidth(1)
        trex_mc_tot.SetFillStyle(0)
        trex_mc_tot.SetLineColor(ROOT.kWhite)
        trex_mc_stat_err = deepcopy(file_temp.Get("g_totErr_postFit"))
        trex_mc_stat_err_only = ROOT.TGraphAsymmErrors()
        for i in range(trex_mc_stat_err.GetN()):
            y = trex_mc_stat_err.GetY()[i]
            x = trex_mc_tot.GetBinCenter(i + 1)
            if y < 1e-10:
                y = 0
                yerr_low = 0
                yerr_high = 0
                trex_mc_stat_err.GetEYlow()[i] = yerr_low
                trex_mc_stat_err.GetEYhigh()[i] = yerr_high
                trex_mc_tot.SetBinContent(i + 1, y)
            xerr_low = trex_mc_tot.GetBinWidth(i + 1) / 2.
            xerr_high = trex_mc_tot.GetBinWidth(i + 1) / 2.
            yerr_low = trex_mc_stat_err.GetEYlow()[i]
            yerr_high = trex_mc_stat_err.GetEYhigh()[i]
            trex_mc_stat_err_only.SetPoint(trex_mc_stat_err_only.GetN(), x, 1)
            if y > 0:
                trex_mc_stat_err_only.SetPointError(trex_mc_stat_err_only.GetN() - 1, xerr_low, xerr_high, yerr_low / y, yerr_high / y)
            else:
                trex_mc_stat_err_only.SetPointError(trex_mc_stat_err_only.GetN() - 1, xerr_low, xerr_high, 0, 0)
        trex_mc_stat_err.SetFillColor(ROOT.kGray + 3)
        trex_mc_stat_err.SetFillStyle(3354)
        trex_mc_stat_err_only.SetFillColor(ROOT.kGray + 3)
        trex_mc_stat_err_only.SetFillStyle(3354)
        for s in mc_map:
            name = s.shortName
            if s.shortName == "MockMC":
                name += "_0tag" if "0tag" in channel.name else "_1tag"
            # if s.shortName == "Wjets_emu_Matched":
            #     name += "_OS" if "OS" in channel.name else "_SS"
            h_temp = deepcopy(file_temp.Get(f"h_{name}_postFit"))
            for i in range(1, h_temp.GetNbinsX() + 1):
                if h_temp.GetBinContent(i) < 1e-10:
                    h_temp.SetBinContent(i, 0)
            mc_map[s] = h_temp
    file_temp.Close()
    return h_data_trex, trex_mc_tot, trex_mc_stat_err, trex_mc_stat_err_only


def save_to_trex_file(trex_folder: str, channel: channel.Channel, var: variable.Variable,
                      h_data: ROOT.TH1, mc_map: MC_Map, trex_histograms: Dict, sys: str = None, affecting: str = ""):
    logging.info(f"Saving histograms to trex root files for channel {channel} for sys {sys}")
    out_name = f"{channel.name}_{var.name}"
    if sys:
        out_name += f"_{sys}"
    if not sys:
        data_file = ROOT.TFile(trex_histograms["Data"], "UPDATE")
        data_file.cd()
        h_data.Write(out_name)
        data_file.Close()
    for s in mc_map:
        if affecting and s.shortName not in affecting:
            continue
        if s.makeGhostSample:
            if not sys:
                h_ghost = mc_map[s].Clone(f"{mc_map[s].GetName()}_GHOST")
                # minimum = h_ghost.GetMinimum()
                for i in range(0, h_ghost.GetNbinsX() + 2):
                    h_ghost.SetBinContent(i, 10000.)
                    h_ghost.SetBinError(i, 0)
                out_file = ROOT.TFile(trex_histograms[s.shortName + "_GHOST"], "UPDATE")
                h_ghost.Write(out_name)
                out_file.Close()
            out_file = ROOT.TFile(trex_histograms[s.shortName], "UPDATE")
            out_file.cd()
            h_offset = mc_map[s].Clone(f"{mc_map[s].GetName()}_OFFSET")
            # for i in range(0, h_offset.GetNbinsX() + 2):
            #     h_offset.SetBinContent(i, h_offset.GetBinContent(i) + 10000.)
            h_offset.Write(out_name)
            out_file.Close()
        else:
            out_file = ROOT.TFile(trex_histograms[s.shortName], "UPDATE")
            out_file.cd()
            mc_map[s].Write(out_name)
            out_file.Close()


def save_histograms_to_trex_file(trex_folder: str, channel: channel.Channel, var: variable.Variable,
                                 h: ROOT.TH1, s: sample.Sample, trex_histograms: Dict, sys: str = None):
    logging.info(f"Saving histogram {h} to trex root files for channel {channel} for sys {sys}")
    out_name = f"{channel.name}_{var.name}"
    if sys:
        out_name += f"_{sys}"
    out_file = ROOT.TFile(trex_histograms[s.shortName], "UPDATE")
    out_file.cd()
    h.Write(out_name)
    out_file.Close()


def save_to_file(out_file_name: str, channel: channel.Channel, var: variable.Variable, h_data: ROOT.TH1, mc_map: MC_Map):
    logging.info(f"Saving histograms to root file for channel {channel}")
    out_file = ROOT.TFile(out_file_name, "UPDATE")
    out_file.cd()
    if h_data:
        h_data.Write()
    for s in mc_map:
        out_name = mc_map[s].GetName()
        if " | " in out_name:
            out_name_split = out_name.split(" | ")
            out_name = f"{out_name_split[0]}_{channel.name}_{var.name}"
        mc_map[s].Write(out_name)
    out_file.Close()


def read_samples(conf: globalConfig.GlobalConfig, reader: inputDataReader.InputDataReader,
                 c: channel.Channel, v: variable.Variable, samples: list,
                 fit: likelihoodFit.LikelihoodFit = None, force_positive: bool = False,
                 sys: str = None, affecting: list = None, fallback: MC_Map = None) -> MC_Map:
    mc_map = {}
    for s in samples:
        if sys and "MockMC" in s.shortName:
            h_nominal = fallback[s]
            h = h_nominal.Clone(f"{h_nominal.GetName()}_{sys}")
            mc_map[s] = h
            continue
        if sys and fallback and affecting:
            if s.shortName not in affecting:
                logging.info(f"Fallback to nominal histogram for {s.shortName} and sys {sys}")
                if s not in fallback:
                    logging.info(f"Fallback histogram not found, continuing...")
                    continue
                h_nominal = fallback[s]
                h = h_nominal.Clone(f"{h_nominal.GetName()}_{sys}")
                mc_map[s] = h
                continue

        # read MC histogram
        h = reader.get_histogram(s, c, v, force_positive, sys)
        if not h:
            continue
        if sys:
            h = h.Clone(f"{h.GetName()}_{sys}")
        mc_map[s] = h

        # scale histogram if performed likelihood fit
        if fit:
            h.Scale(fit.result[s.shortName][0])
    return mc_map


def get_lumi(lumi_string):
    lumi = 0
    if "2015" in lumi_string:
        lumi += 3219.56
    if "2016A-B" in lumi_string:
        lumi += 2289.2
    if "2016C-L" in lumi_string:
        lumi += 30698.9
    if "2016" in lumi_string and "A-B" not in lumi_string and "C-L" not in lumi_string:
        lumi += 32988.1
    if "2017" in lumi_string:
        lumi += 44307.4
    if "2018" in lumi_string:
        lumi += 58450.1
    return lumi


def get_variables(options, conf, reader, channel, sample=None):
    variables_conf = conf.get_variable_names()
    variables = []
    if options.vars:
        variables = options.vars.split(",")
    else:
        if sample:
            test_sample = sample
        else:
            test_sample = conf.data
        if channel.make_plots:
            for v in reader.find_variables(test_sample, channel):
                variables += [v]

    out = []
    for v in variables_conf:
        if v in variables:
            out += [v]

    return out


def get_maximum(h, x1, x2):
    out = h.GetBinContent(h.FindBin(x1))
    for i in range(1, h.GetNbinsX() + 1):
        if h.GetBinCenter(i) < x1 or h.GetBinCenter(i) > x2:
            continue
        if h.GetBinContent(i) > out:
            out = h.GetBinContent(i)
    return out


def set_to_positive(h, sys=""):
    for i in range(0, h.GetNbinsX() + 2):
        if h.GetBinContent(i) <= 0:
            if not sys:
                h.SetBinContent(i, 1e-5)
            else:
                h.SetBinContent(i, 1e-6)
            # h.SetBinError(i, 1e-5)


def fit_histogram(h, fit):
    func = ROOT.TF1(f"{h.GetName()}_fit", fit)
    h.Fit(func, "0")
    print(func)
    for i in range(1, h.GetNbinsX() + 1):
        h.SetBinContent(i, func.Integral(h.GetBinLowEdge(i), h.GetBinLowEdge(i) + h.GetBinWidth(i)) / h.GetBinWidth(i))
        h.SetBinError(i, 0)
    logging.info(f"Fit Chi2 prob {ROOT.TMath.Prob(func.GetChisquare(), func.GetNDF())}")
    return h


def set_errors_to_zero(h):
    for i in range(0, h.GetNbinsX() + 2):
        h.SetBinError(i, 0.)


def set_under_over_flow(h: ROOT.TH1, x_range: list, do_overflow: bool, do_underflow: bool):
    N = h.GetNbinsX()
    width = h.GetBinWidth(1)
    if not x_range:
        x_range = [h.GetBinCenter(1) - h.GetBinWidth(1) * 0.5, h.GetBinCenter(N) + h.GetBinWidth(N) * 0.5]
    x_range_bins = [h.FindBin(x_range[0] + width * 0.1), h.FindBin(x_range[1] - width * 0.1)]

    h_new = ROOT.TH1F(f"{h.GetName()}_rebined", f"{h.GetName()}_rebined", x_range_bins[1] - x_range_bins[0] + 1, x_range[0], x_range[1])

    err0 = c_double()
    err1 = c_double()
    errN = c_double()
    errN1 = c_double()

    val0 = h.IntegralAndError(0, x_range_bins[0] - 1, err0)
    val1 = h.IntegralAndError(x_range_bins[0], x_range_bins[0], err1)
    valN = h.IntegralAndError(x_range_bins[1], x_range_bins[1], errN)
    valN1 = h.IntegralAndError(x_range_bins[1] + 1, N + 1, errN1)

    for i in range(0, x_range_bins[0]):
        h.SetBinContent(i, 0)
        h.SetBinError(i, 0)
    for i in range(x_range_bins[1] + 1, N + 2):
        h.SetBinContent(i, 0)
        h.SetBinError(i, 0)

    if do_underflow:
        h.SetBinContent(x_range_bins[0], val0 + val1)
        h.SetBinError(x_range_bins[0], (err0.value**2 + err1.value**2)**(0.5))
    if do_overflow:
        h.SetBinContent(x_range_bins[1], valN + valN1)
        h.SetBinError(x_range_bins[1], (errN.value**2 + errN1.value**2)**(0.5))

    j = 1
    for i in range(x_range_bins[0], x_range_bins[1] + 1):
        h_new.SetBinContent(j, h.GetBinContent(i))
        h_new.SetBinError(j, h.GetBinError(i))
        j += 1

    h_new.SetName(h.GetName())
    h_new.SetTitle(h.GetTitle())

    return h_new


def rebin_histogram(h: ROOT.TH1, v: variable.Variable, extra_rebin: int = 1):
    rebin = v.rebin
    if rebin and v.allow_rebin and not v.dstar_tail_rebin:
        h.Rebin(int(rebin * extra_rebin))
        if v.allow_rebin:
            return set_under_over_flow(h, v.x_range, v.do_overflow, v.do_underflow)
        else:
            return h
    elif rebin and v.allow_rebin and v.dstar_tail_rebin and v.name == "Dmeson_mdiff":
        if extra_rebin < 1.5:
            x_var_list = [135., 140., 141., 142., 143., 144., 145., 146., 147., 148., 149., 150., 153., 156., 159., 162., 165., 168., 171., 174., 177., 180.]
            x_var = array.array('d', x_var_list)
            h_out = h.Rebin(len(x_var_list) - 1, "h_out", x_var)
            return h_out
        elif extra_rebin > 1.5:
            x_var_list = [140., 142., 144., 146., 148., 150., 156., 162., 168., 174., 180.]
            x_var = array.array('d', x_var_list)
            h_out = h.Rebin(len(x_var_list) - 1, "h_out", x_var)
            return h_out
    else:
        return h


def make_stat_err(h: ROOT.TH1) -> List[Union[ROOT.TGraphErrors, ROOT.TGraphErrors]]:
    gr = ROOT.TGraphAsymmErrors()
    gr_err_only = ROOT.TGraphAsymmErrors()
    for i in range(0, h.GetNbinsX() + 2):
        x = h.GetBinCenter(i)
        y = h.GetBinContent(i)
        gr.SetPoint(gr.GetN(), x, y)
        gr.SetPointError(gr.GetN() - 1, (h.GetBinWidth(i) / 2.), (h.GetBinWidth(i) / 2.), h.GetBinError(i), h.GetBinError(i))
        gr_err_only.SetPoint(gr_err_only.GetN(), x, 1)
        if y:
            gr_err_only.SetPointError(gr_err_only.GetN() - 1, (h.GetBinWidth(i) / 2.), (h.GetBinWidth(i) / 2.), h.GetBinError(i) / y, h.GetBinError(i) / y)
        else:
            gr_err_only.SetPointError(gr_err_only.GetN() - 1, (h.GetBinWidth(i) / 2.), (h.GetBinWidth(i) / 2.), 0, 0)
    gr.SetFillColor(ROOT.kGray + 3)
    gr.SetFillStyle(3354)
    gr_err_only.SetFillColor(ROOT.kGray + 3)
    gr_err_only.SetFillStyle(3354)
    return gr, gr_err_only


def make_stat_err_and_nominal(h: ROOT.TH1) -> List[Union[ROOT.TGraphErrors, ROOT.TGraphErrors, ROOT.TGraphErrors]]:
    gr_val_only = ROOT.TGraphErrors()
    gr = ROOT.TGraphAsymmErrors()
    gr_err_only = ROOT.TGraphAsymmErrors()
    for i in range(0, h.GetNbinsX() + 2):
        x = h.GetBinCenter(i)
        y = h.GetBinContent(i)
        gr_val_only.SetPoint(gr_val_only.GetN(), x, y)
        gr_val_only.SetPointError(gr_val_only.GetN() - 1, (h.GetBinWidth(i) / 2.), 0)
        gr.SetPoint(gr.GetN(), x, y)
        gr.SetPointError(gr.GetN() - 1, (h.GetBinWidth(i) / 2.), (h.GetBinWidth(i) / 2.), h.GetBinError(i), h.GetBinError(i))
        gr_err_only.SetPoint(gr_err_only.GetN(), x, 1)
        if y:
            gr_err_only.SetPointError(gr_err_only.GetN() - 1, (h.GetBinWidth(i) / 2.), (h.GetBinWidth(i) / 2.), h.GetBinError(i) / y, h.GetBinError(i) / y)
        else:
            gr_err_only.SetPointError(gr_err_only.GetN() - 1, (h.GetBinWidth(i) / 2.), (h.GetBinWidth(i) / 2.), 0, 0)
    gr.SetFillColor(ROOT.kGray + 3)
    gr.SetFillStyle(3354)
    gr_err_only.SetFillColor(ROOT.kGray + 3)
    gr_err_only.SetFillStyle(3354)
    return gr_val_only, gr, gr_err_only


def make_sys_err(h: ROOT.TH1, h_sys: List) -> List[Union[ROOT.TGraphErrors, ROOT.TGraphErrors]]:
    gr = ROOT.TGraphAsymmErrors()
    gr_err_only = ROOT.TGraphAsymmErrors()
    for i in range(0, h.GetNbinsX() + 2):
        x = h.GetBinCenter(i)
        y = h.GetBinContent(i)
        sys = [h_sys.GetBinContent(i) for h_sys in h_sys]
        y_err_up = 0
        y_err_dn = 0
        for y_sys in sys:
            if y_sys >= y:
                err = y_sys - y
                y_err_up += err**2
            else:
                err = y - y_sys
                y_err_dn += err**2
        y_err_up = (y_err_up)**(0.5)
        y_err_dn = (y_err_dn)**(0.5)
        gr.SetPoint(gr.GetN(), x, y)
        gr.SetPointError(gr.GetN() - 1, (h.GetBinWidth(i) / 2.), (h.GetBinWidth(i) / 2.), y_err_dn, y_err_up)
        gr_err_only.SetPoint(gr_err_only.GetN(), x, 1)
        if y:
            gr_err_only.SetPointError(gr_err_only.GetN() - 1, (h.GetBinWidth(i) / 2.), (h.GetBinWidth(i) / 2.), y_err_dn / y, y_err_up / y)
        else:
            gr_err_only.SetPointError(gr_err_only.GetN() - 1, (h.GetBinWidth(i) / 2.), (h.GetBinWidth(i) / 2.), 0, 0)
    gr.SetFillColor(ROOT.kGray + 3)
    gr.SetFillStyle(3354)
    gr_err_only.SetFillColor(ROOT.kGray + 3)
    gr_err_only.SetFillStyle(3354)
    return gr, gr_err_only


def make_empty_error_bands(h: ROOT.TH1) -> List[Union[ROOT.TGraphErrors, ROOT.TGraphErrors]]:
    gr = ROOT.TGraphAsymmErrors()
    gr_err_only = ROOT.TGraphAsymmErrors()
    for i in range(0, h.GetNbinsX() + 2):
        x = h.GetBinCenter(i)
        y = h.GetBinContent(i)
        gr.SetPoint(gr.GetN(), x, y)
        gr.SetPointError(gr.GetN() - 1, (h.GetBinWidth(i) / 2.), (h.GetBinWidth(i) / 2.), 0, 0)
        gr_err_only.SetPoint(gr_err_only.GetN(), x, 1)
        gr_err_only.SetPointError(gr_err_only.GetN() - 1, (h.GetBinWidth(i) / 2.), (h.GetBinWidth(i) / 2.), 0, 0)
    return gr, gr_err_only


def make_pdf_err(h: ROOT.TH1, h_var: List, pdfstring: str, norm: bool = False) -> List[Union[ROOT.TGraphErrors, ROOT.TGraphErrors]]:
    # h_var is expected to have only the PDF variations and no QCD parameter variations!!
    if not len(h_var):
        return make_empty_error_bands(h)

    pdfset = lhapdf.getPDFSet(pdfstring)
    pdfset.mkPDFs()
    nBins = h.GetNbinsX()
    if norm:
        nom_sum = h_var[0].GetSum()
        for hist in h_var:
            hist.Scale(nom_sum / hist.GetSum())

    xval = ROOT.std.vector("float")(nBins)
    yval = ROOT.std.vector("float")(nBins)
    exh = ROOT.std.vector("float")(nBins)
    exl = ROOT.std.vector("float")(nBins)
    for a in range(nBins):
        xval[a] = h.GetBinCenter(a + 1)
        yval[a] = h.GetBinContent(a + 1)
        exh[a] = h.GetBinWidth(a + 1) / 2.
        exl[a] = h.GetBinWidth(a + 1) / 2.

    exh = np.asarray(exh)
    exl = np.asarray(exl)
    xval = np.asarray(xval)
    yval = np.asarray(yval)
    pdfvars = np.zeros((len(h_var), nBins))
    for i in range(len(h_var)):
        pdfvars[i, :] = np.array([h_var[i].GetBinContent(binNo) for binNo in range(1, nBins + 1)])
    pdferrplus = np.array([pdfset.uncertainty(pdfvars[:, i]).errplus for i in range(nBins)])
    pdferrminus = np.array([pdfset.uncertainty(pdfvars[:, i]).errminus for i in range(nBins)])
    exh = np.hstack((0, exh, 0))
    exl = np.hstack((0, exl, 0))
    xval = np.hstack((xval[0], xval, xval[-1]))
    yval = np.hstack((yval[0], yval, yval[-1]))
    pdferrplus = np.hstack(([0], pdferrplus, [0]))
    pdferrminus = np.hstack(([0], pdferrminus, [0]))
    with np.errstate(divide='ignore', invalid='ignore'):
        pdfeyh_o = pdferrplus / yval
        pdfeyl_o = pdferrminus / yval
        pdfeyh_o[yval == 0.] = 0.
        pdfeyl_o[yval == 0.] = 0.
    pdfeyh = array.array('d', pdferrplus)
    pdfeyl = array.array('d', pdferrminus)
    pdfeyh_o = array.array('d', pdfeyh_o)
    pdfeyl_o = array.array('d', pdfeyl_o)
    exh = array.array('d', exh)
    exl = array.array('d', exl)
    unity = array.array('d', np.ones((nBins + 2,)))
    xval = array.array('d', xval)
    yval = array.array('d', yval)
    gr = ROOT.TGraphAsymmErrors(nBins + 2, xval, yval, exl, exh, pdfeyl, pdfeyh)
    gr_err_only = ROOT.TGraphAsymmErrors(nBins + 2, xval, unity, exl, exh, pdfeyl_o, pdfeyh_o)
    gr.SetFillColor(ROOT.kGray + 3)
    gr.SetFillStyle(3354)
    gr_err_only.SetFillColor(ROOT.kGray + 3)
    gr_err_only.SetFillStyle(3354)
    return gr, gr_err_only


def make_minmax_err(h: ROOT.TH1, h_var: List, Norm: bool = False) -> List[Union[ROOT.TGraphErrors, ROOT.TGraphErrors]]:
    # h_var is expected to have only the QCD parameter variations!!
    if not len(h_var):
        return make_empty_error_bands(h)

    nBins = h.GetNbinsX()

    if Norm:
        nom_sum = h_var[0].GetSum()
        for hist in h_var:
            hist.Scale(nom_sum / hist.GetSum())

    xval = ROOT.std.vector("float")(nBins)
    yval = ROOT.std.vector("float")(nBins)
    exh = ROOT.std.vector("float")(nBins)
    exl = ROOT.std.vector("float")(nBins)

    eyh = eyl = np.zeros((nBins,))

    for a in range(nBins):
        xval[a] = h.GetBinCenter(a + 1)
        yval[a] = h.GetBinContent(a + 1)
        exh[a] = h.GetBinWidth(a + 1) / 2.
        exl[a] = h.GetBinWidth(a + 1) / 2.

    exh = np.asarray(exh)
    exl = np.asarray(exl)

    xval = np.asarray(xval)
    yval = np.asarray(yval)

    for i in range(1, h.GetNbinsX() + 1):
        bin_list = []
        for hist in h_var:
            bin_list.append(hist.GetBinContent(i))
        if len(bin_list) != 0:
            eyh[i - 1] = np.abs(np.max(bin_list) - yval[i - 1])
            eyl[i - 1] = np.abs(np.min(bin_list) - yval[i - 1])
        else:
            eyh[i - 1] = 0
            eyl[i - 1] = 0

    exh = np.hstack((0, exh, 0))
    exl = np.hstack((0, exl, 0))
    xval = np.hstack((xval[0], xval, xval[-1]))
    yval = np.hstack((yval[0], yval, yval[-1]))
    eyh = np.hstack((0, eyh, 0))
    eyl = np.hstack((0, eyl, 0))

    with np.errstate(divide='ignore', invalid='ignore'):
        eyh_o = eyh / yval
        eyl_o = eyl / yval
        eyh_o[yval == 0.] = 0.
        eyl_o[yval == 0.] = 0.

    eyh = array.array('d', eyh)
    eyl = array.array('d', eyl)

    eyh_o = array.array('d', eyh_o)
    eyl_o = array.array('d', eyl_o)

    exh = array.array('d', exh)
    exl = array.array('d', exl)

    xval = array.array('d', xval)
    yval = array.array('d', yval)
    unity = array.array('d', np.ones((nBins + 2,)))
    gr = ROOT.TGraphAsymmErrors(nBins + 2, xval, yval, exl, exl, eyl, eyh)
    gr_err_only = ROOT.TGraphAsymmErrors(nBins + 2, xval, unity, exl, exl, eyl_o, eyh_o)

    gr.SetFillColor(ROOT.kGray + 3)
    gr.SetFillStyle(3354)
    gr_err_only.SetFillColor(ROOT.kGray + 3)
    gr_err_only.SetFillStyle(3354)
    return gr, gr_err_only


def combine_error(gr1: ROOT.TGraphAsymmErrors, gr2: ROOT.TGraphAsymmErrors) -> ROOT.TGraphAsymmErrors:
    gr = ROOT.TGraphAsymmErrors()
    for i in range(0, gr1.GetN()):
        x = gr1.GetX()[i]
        y = gr1.GetY()[i]
        err_x_dn = gr1.GetEXlow()[i]
        err_x_up = gr1.GetEXhigh()[i]
        err_y_1_dn = gr1.GetEYlow()[i]
        err_y_1_up = gr1.GetEYhigh()[i]
        err_y_2_dn = gr2.GetEYlow()[i]
        err_y_2_up = gr2.GetEYhigh()[i]
        gr.SetPoint(i, x, y)
        gr.SetPointError(i, err_x_dn, err_x_up, (err_y_1_dn**2 + err_y_2_dn**2)**(0.5), (err_y_1_up**2 + err_y_2_up**2)**(0.5))
    gr.SetFillColor(ROOT.kBlue - 4)
    gr.SetFillStyle(3345)
    return gr


def combine_error_multiple(asym_list):

    gr1 = asym_list[0]
    gr = ROOT.TGraphAsymmErrors()
    for i in range(0, gr1.GetN()):
        x = gr1.GetX()[i]
        y = gr1.GetY()[i]
        err_x_dn = gr1.GetEXlow()[i]
        err_x_up = gr1.GetEXhigh()[i]
        err_y_dn = np.sum(np.array([asym.GetEYlow()[i] for asym in asym_list])**2)**0.5
        err_y_up = np.sum(np.array([asym.GetEYhigh()[i] for asym in asym_list])**2)**0.5
        gr.SetPoint(i, x, y)
        gr.SetPointError(i, err_x_dn, err_x_up, err_y_dn, err_y_up)

    gr.SetFillColor(ROOT.kGray + 3)
    gr.SetFillStyle(3354)
    # gr.SetFillColor(ROOT.kBlue - 4)
    # gr.SetFillStyle(3345)
    return gr


def make_ratio(data: ROOT.TH1, mc_tot: ROOT.TH1) -> ROOT.TH1:
    h_ratio = data.Clone(f"{data.GetName()}_over_{mc_tot.GetName()}")
    for i in range(0, h_ratio.GetNbinsX() + 2):
        y_mc = mc_tot.GetBinContent(i)
        y_data = data.GetBinContent(i)
        err_data = data.GetBinError(i)
        if y_mc and y_data:
            ratio = data.GetBinContent(i) / mc_tot.GetBinContent(i)
            h_ratio.SetBinContent(i, ratio)
            h_ratio.SetBinError(i, ratio * (err_data / y_data))
        else:
            h_ratio.SetBinContent(i, 0)
            h_ratio.SetBinError(i, 0)
    return h_ratio


def make_stack(samples: List, mc_map: MC_Map):
    hs = ROOT.THStack()
    for s in reversed(samples):
        if s not in mc_map.keys():
            continue
        h = mc_map[s]
        logger.debug(f"adding {h}")
        hs.Add(h)
    return hs


def sys_graph_to_hists(gr: ROOT.TGraphErrors, nominal: ROOT.TH1, name: str) -> List[Union[ROOT.TH1, ROOT.TH1]]:
    h_up = nominal.Clone(f"{nominal.GetName()}_{name}_up")
    h_dn = nominal.Clone(f"{nominal.GetName()}_{name}_dn")
    for i in range(0, nominal.GetNbinsX() + 2):
        y = nominal.GetBinContent(i)
        h_up.SetBinContent(i, y + gr.GetEYhigh()[i])
        h_dn.SetBinContent(i, y - gr.GetEYlow()[i])
    return h_up, h_dn


def make_mc_tot(hs: ROOT.THStack, name: str) -> ROOT.TH1:
    logger.debug("In make mc tot")
    h = hs.GetStack().Last().Clone(name)
    h.SetLineWidth(1)
    h.SetFillStyle(0)
    h.SetLineColor(ROOT.kWhite)
    return h


def get_fraction_histogram(h1, h2):
    h = h1.Clone(f"{h1.GetName()}_fraction")
    h.Divide(h2)
    h.SetFillStyle(0)
    h.SetLineWidth(2)
    h.SetLineColor(ROOT.kBlack)
    gr_err, _ = make_stat_err(h)
    return h, gr_err


def get_samples(conf: globalConfig.GlobalConfig, channel: channel.Channel) -> List:
    if channel.samples:
        samples = [conf.get_sample(s) for s in channel.samples]
    else:
        samples = conf.get_mc()
    return samples


def mass_fit(conf: globalConfig.GlobalConfig, h: ROOT.TH1, channel: channel.Channel, label: str) -> massFit.MassFit:
    if channel.mass_fit:
        variable = conf.get_var(channel.mass_fit["var"])
        fit = massFit.MassFit(channel, variable, h, label, os.path.join(conf.config_name, f"{channel.name}_mass_fit"))
        fit.fit()


def likelihood_fit(conf: globalConfig.GlobalConfig, reader: inputDataReader.InputDataReader,
                   channel: channel.Channel, samples: List) -> likelihoodFit.LikelihoodFit:
    if channel.likelihood_fit:
        # perform the fit
        logger.debug(f"performing likelihood fit with configuration {channel.likelihood_fit}")
        fit_var = conf.get_var(channel.likelihood_fit['variable'])
        fit_range = channel.likelihood_fit['range']
        h_data = reader.get_histogram(conf.get_data(), channel, fit_var)
        mc_map = {}
        for s in samples:
            h = reader.get_histogram(s, channel, fit_var)
            mc_map[s] = h
        fit = likelihoodFit.LikelihoodFit(channel, h_data, mc_map, fit_range, conf.config_name)
        fit.fit()

        # extrapolate to SR
        if 'extrapolated_region' in channel.likelihood_fit:
            logger.debug(f"extrapolating fit to signal region {channel.likelihood_fit['extrapolated_region']}")
            channel_SR = conf.get_channel(channel.likelihood_fit['extrapolated_region'])
            samples_SR = [conf.get_sample(s) for s in channel_SR.samples]
            h_data_SR = reader.get_histogram(conf.get_data(), channel_SR, fit_var)
            mc_map_SR = {}
            for s in samples_SR:
                h = reader.get_histogram(s, channel_SR, fit_var)
                mc_map_SR[s] = h
            for s in mc_map_SR:
                if 'Multijet' in s.name:
                    continue
                h = mc_map_SR[s]
                h.Scale(fit.result[s.shortName][0])
                print(s.name, " ", s.shortName, " ", fit.result[s.shortName][0])
            sf_qcd_SR = scale_multijet_histogram(h_data_SR, mc_map_SR, fit_range)
            # if sf_qcd_SR < 0:
            #     sf_qcd_SR = 0
            for s in mc_map_SR:
                if 'Multijet' in s.name:
                    fit.result[s.name] = [sf_qcd_SR, 0]
                    break

        # return fit
        fit.save_results()
        return fit
    else:
        return None


def scale_multijet_histogram(data: ROOT.TH1, mc_map: MC_Map, fit_range: list):
    bin1 = data.FindBin(fit_range[0])
    bin2 = data.FindBin(fit_range[1])
    integral_data = data.Integral(bin1, bin2)
    integral_ewk = 0
    integral_qcd = 0
    for s in mc_map:
        h = mc_map[s]
        if 'Multijet' not in s.name:
            integral_ewk += h.Integral(bin1, bin2)
        else:
            integral_qcd = h.Integral(bin1, bin2)
    return (integral_data - integral_ewk) / integral_qcd


def average_content(h, i):
    val = 0
    N = 0
    if i > 1:
        val += h.GetBinContent(i - 1)
        N += 1
    if i < h.GetNbinsX():
        val += h.GetBinContent(i + 1)
        N += 1
    return val / N


def replace_sample(conf: globalConfig.GlobalConfig, mc_map: MC_Map, reader: inputDataReader.InputDataReader,
                   c: channel.Channel, var: variable.Variable, sample: str, channel: str, mc_map_sys: Dict[str, MC_Map] = None):
    channel_replacement = conf.get_channel(channel)
    assert len(channel_replacement.samples) == 1
    sample_replacement = conf.get_sample(channel_replacement.samples[0])
    sample_current = [conf.get_sample(s) for s in c.samples if conf.get_sample(s).shortName == sample]
    print(sample_current)
    assert len(sample_current) == 1
    sample_current = sample_current[0]
    h_replacement = reader.get_histogram(sample_replacement, channel_replacement, var, channel_replacement.force_positive)
    h_current = mc_map[sample_current]
    h_current = h_current.Clone(f"{h_current.GetName()}_temp_clone")
    h_replacement.Scale(h_current.Integral(0, h_current.GetNbinsX() + 1) / h_replacement.Integral(0, h_replacement.GetNbinsX() + 1))
    if h_replacement.GetNbinsX() > h_current.GetNbinsX():
        logging.warning(f"Replacement sample has more bins. Rebinning {h_replacement.GetNbinsX()} -> {h_current.GetNbinsX()}")
        h_replacement.Rebin(int(h_replacement.GetNbinsX() / h_current.GetNbinsX()))
    mc_map[sample_current] = h_replacement
    if mc_map_sys:
        # smooth any bins with very low bin content in the nominal histogram
        for i in range(0, h_current.GetNbinsX() + 2):
            if abs(h_current.GetBinContent(i)) < 1e-3:
                h_current.SetBinContent(i, average_content(h_current, i))
        for group, systematics in mc_map_sys.items():
            for syst, map_sys in systematics.items():
                h_sys_replaced = map_sys[sample_current].Clone(f"{map_sys[sample_current].GetName()}_replaced")
                h_sys_replaced.Add(h_current, -1)

                # use relative uncertaitny instead of absolute in these cases
                if "PROXY_NORM" in syst:
                    h_sys_replaced.Divide(h_current)
                    h_sys_replaced.Multiply(h_replacement)

                h_sys_replaced.Add(h_replacement)
                map_sys[sample_current] = h_sys_replaced


def make_canvas(h: ROOT.TH1, v: variable.Variable, c: channel.Channel,
                x: float = 800., y: float = 600., y_split: float = 0.30,
                fit: likelihoodFit.LikelihoodFit = None, scale_factors: dict = None,
                events: str = "Entries", sys: str = None) -> ROOT.TCanvas:
    canv = canvas.Canvas2(c, v, x, y, y_split, fit, scale_factors, sys=sys)
    canv.construct(h, events=events)
    return canv


def make_canvas_unfold(h: ROOT.TH1, v: variable.Variable, c: channel.Channel,
                       x: float = 800., y: float = 600., y_split: float = 0.30,
                       fit: likelihoodFit.LikelihoodFit = None, scale_factors: dict = None,
                       events: str = "Entries", sys: str = None) -> ROOT.TCanvas:
    canv = canvas.CanvasCrossSection(c, v, x, y, y_split, fit, scale_factors, sys=sys)
    canv.construct(h, events=events)
    return canv


def make_canvas_mc_ratio(h: ROOT.TH1, v: variable.Variable, c: channel.Channel, ratio_title: str,
                         x: float = 800., y: float = 600., y_split: float = 0.30,
                         ratio_range: list = [0.01, 1.99], events: str = "Entries") -> ROOT.TCanvas:
    canv = canvas.CanvasMCRatio(c, v, ratio_title, x, y, y_split, ratio_range)
    canv.construct(h, events=events)
    return canv
