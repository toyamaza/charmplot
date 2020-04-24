from charmplot.common import canvas
from charmplot.common import likelihoodFit
from charmplot.common import massFit
from charmplot.control import channel
from charmplot.control import globalConfig
from charmplot.control import inputDataReader
from charmplot.control import sample
from charmplot.control import variable
from typing import Dict, List, Union
import json
import logging
import os
import ROOT
import sys

MC_Map = Dict[sample.Sample, ROOT.TH1]

logger = logging.getLogger(__name__)


def get_mc_min(mc_map: MC_Map, samples: list):
    for s in reversed(samples):
        if s not in mc_map.keys():
            continue
        return mc_map[s]


def save_to_file(out_file_name: str, channel: channel.Channel, var: variable.Variable, h_data: ROOT.TH1, mc_map: MC_Map):
    logging.info(f"Saving histograms to root file for channel {channel}")
    out_file = ROOT.TFile(out_file_name, "UPDATE")
    out_file.cd()
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
                 fit: likelihoodFit.LikelihoodFit=None) -> MC_Map:
    mc_map = {}
    for s in samples:
        # read MC histogram
        h = reader.get_histogram(s, c, v)
        if not h:
            continue
        mc_map[s] = h

        # scale histogram if performed likelihood fit
        if fit:
            h.Scale(fit.result[s.shortName][0])
    return mc_map


def get_lumi(lumi_string):
    lumi = 0
    if "2015" in lumi_string:
        lumi += 3219.56
    if "2016" in lumi_string:
        lumi += 32988.1
    if "2017" in lumi_string:
        lumi += 43587.3
    if "2018" in lumi_string:
        lumi += 58450.1
    return lumi


def get_variables(options, conf, reader, channel, sample=None):
    variables_conf = conf.get_variable_names()
    variables = []
    picked_variables = []
    if options.vars:
        picked_variables = options.vars.split(",")
    if sample:
        test_sample = sample
    else:
        test_sample = conf.data
    for v in reader.find_variables(test_sample, channel):
        if picked_variables and v not in picked_variables:
            continue
        variables += [v]

    out = []
    for v in variables_conf:
        if v in variables:
            out += [v]
            variables.remove(v)

    return out + variables


def get_maximum(h, x1, x2):
    out = h.GetBinContent(h.FindBin(x1))
    for i in range(1, h.GetNbinsX() + 1):
        if h.GetBinCenter(i) < x1 or h.GetBinCenter(i) > x2:
            continue
        if h.GetBinContent(i) > out:
            out = h.GetBinContent(i)
    return out


def set_to_positive(h):
    for i in range(0, h.GetNbinsX() + 2):
        if h.GetBinContent(i) < 0:
            h.SetBinContent(i, 0)


def set_under_over_flow(h: ROOT.TH1, x_range: list):
    N = h.GetNbinsX()
    if not x_range:
        x_range = [h.GetBinCenter(1) - h.GetBinWidth(1) * 0.5, h.GetBinCenter(N) + h.GetBinWidth(N) * 0.5]
    x_range_bins = [h.FindBin(x_range[0]), h.FindBin(x_range[1] - h.GetBinWidth(N) * 0.1)]

    val0 = ROOT.Double()
    err0 = ROOT.Double()
    val1 = ROOT.Double()
    err1 = ROOT.Double()
    valN = ROOT.Double()
    errN = ROOT.Double()
    valN1 = ROOT.Double()
    errN1 = ROOT.Double()

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

    h.SetBinContent(x_range_bins[0], val0 + val1)
    h.SetBinError(x_range_bins[0], (err0**2 + err1**2)**(0.5))
    h.SetBinContent(x_range_bins[1], valN + valN1)
    h.SetBinError(x_range_bins[1], (errN**2 + errN1**2)**(0.5))


def rebin_histogram(h: ROOT.TH1, v: variable.Variable, extra_rebin: int = 1):
    rebin = v.rebin
    if rebin:
        h.Rebin(rebin * extra_rebin)
    set_under_over_flow(h, v.x_range)


def make_stat_err(h: ROOT.TH1) -> List[Union[ROOT.TGraphErrors, ROOT.TGraphErrors]]:
    gr = ROOT.TGraphErrors()
    gr_err_only = ROOT.TGraphErrors()
    for i in range(0, h.GetNbinsX() + 2):
        x = h.GetBinCenter(i)
        y = h.GetBinContent(i)
        gr.SetPoint(gr.GetN(), x, y)
        gr.SetPointError(gr.GetN() - 1, (h.GetBinWidth(i) / 2.), h.GetBinError(i))
        gr_err_only.SetPoint(gr_err_only.GetN(), x, 1)
        if y:
            gr_err_only.SetPointError(gr_err_only.GetN() - 1, (h.GetBinWidth(i) / 2.), h.GetBinError(i) / y)
        else:
            gr_err_only.SetPointError(gr_err_only.GetN() - 1, (h.GetBinWidth(i) / 2.), 0)
    gr.SetFillColor(ROOT.kGray + 2)
    gr.SetFillStyle(3354)
    gr_err_only.SetFillColor(ROOT.kGray + 2)
    gr_err_only.SetFillStyle(3354)
    return gr, gr_err_only


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


def make_mc_tot(hs: ROOT.THStack, name: str) -> ROOT.TH1:
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


def read_scale_factors(scale_factor_confing: dict):
    if not scale_factor_confing or 'input_file' not in scale_factor_confing:
        return None
    path = scale_factor_confing['input_file']
    if not os.path.isfile(path):
        logger.warning(f"scale factor file not found: {path}")
        sys.exit(1)
    with open(path, 'r') as json_file:
        logger.debug(f"scale factor file successfully opened: {path}")
        data = json.load(json_file)
        return data


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


def make_canvas(h: ROOT.TH1, v: variable.Variable, c: channel.Channel,
                x: float = 800., y: float = 600., y_split: float = 0.30,
                fit: likelihoodFit.LikelihoodFit = None, scale_factors: dict = None) -> ROOT.TCanvas:
    canv = canvas.Canvas2(c, v, x, y, y_split, fit, scale_factors)
    canv.construct(h)
    return canv


def make_canvas_mc_ratio(h: ROOT.TH1, v: variable.Variable, c: channel.Channel,
                         x: float = 800., y: float = 600., y_split: float = 0.30) -> ROOT.TCanvas:
    canv = canvas.CanvasMCRatio(c, v, x, y, y_split)
    canv.construct(h)
    return canv
