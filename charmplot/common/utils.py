from charmplot.common import canvas
from charmplot.control import channel
from charmplot.control import globalConfig
from charmplot.control import sample
from charmplot.control import variable
from typing import Dict, List, Union
import logging
import ROOT

MC_Map = Dict[sample.Sample, ROOT.TH1]

logger = logging.getLogger(__name__)


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


def get_variables(options, conf, reader, channel):
    variables_conf = conf.get_variable_names()
    variables = []
    picked_variables = []
    if options.vars:
        picked_variables = options.vars.split(",")
    for v in reader.find_variables(conf.data, channel):
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


def set_under_over_flow(h: ROOT.TH1, x_range: list):
    N = h.GetNbinsX()
    if not x_range:
        x_range = [h.GetBinCenter(1) - h.GetBinWidth(1) / 0.5, h.GetBinCenter(N) + h.GetBinWidth(N) / 0.5]
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


def rebin_histogram(h: ROOT.TH1, v: variable.Variable):
    rebin = v.rebin
    if rebin:
        h.Rebin(rebin)
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


def make_stack(conf: globalConfig.GlobalConfig, mc_map: MC_Map):
    hs = ROOT.THStack()
    for s in reversed(conf.get_mc()):
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


def make_canvas(h: ROOT.TH1, v: variable.Variable, c: channel.Channel, r: float = 800 / 800., y_split: float = 0.30) -> List[Union[ROOT.TCanvas, ROOT.TH1]]:
    canv = canvas.Canvas2(c, v, r, y_split)
    canv.construct(h)
    return canv


def make_canvas_mc_ratio(h: ROOT.TH1, v: variable.Variable, c: channel.Channel,
                         r: float = 800 / 800., y_split: float = 0.30) -> List[Union[ROOT.TCanvas, ROOT.TH1]]:
    canv = canvas.CanvasMCRatio(c, v, r, y_split)
    canv.construct(h)
    return canv
