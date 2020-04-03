from charmplot.control import channel
from charmplot.control import sample
from typing import Dict, List
import json
import logging
import os
import ROOT


MC_Map = Dict[sample.Sample, ROOT.TH1]

# logging
logger = logging.getLogger(__name__)


def get_ranged_histogram(h, xmin, xmax):
    first_bin = h.FindBin(xmin)
    last_bin = h.FindBin(xmax)
    if xmax % h.GetBinWidth(last_bin) == 0:
        last_bin -= 1
    nbins = last_bin - first_bin + 1
    if nbins < 1:
        return None
    h_ranged = ROOT.TH1F(f"{h.GetName()}_ranged", "", nbins, 0, nbins)
    h_ranged.Sumw2()
    nbin = 0
    for i in range(first_bin, last_bin + 1):
        nbin += 1
        cont = h.GetBinContent(i)
        if cont <= 0:
            cont = 1e-6
        h_ranged.SetBinContent(nbin, cont)
        h_ranged.SetBinError(nbin, h.GetBinError(i))
    return h_ranged


class LikelihoodFit(object):

    data_integral = 0
    samples = []
    fixed = []
    result = {}
    mc_integral = {}

    def __init__(self, chanenl: channel.Channel, data: ROOT.TH1, mc_map: MC_Map, x_range: List = [], output: str = ""):

        self.chanenl = chanenl
        self.name = chanenl.name
        self.samples = [s for s in mc_map]
        self.output = output

        mc = ROOT.TObjArray(len(self.samples))
        for s in self.samples:
            h = mc_map[s].Clone(f"{mc_map[s].GetName()}_fit")
            if x_range:
                h = get_ranged_histogram(h, x_range[0], x_range[1])
            self.mc_integral[s] = h.GetSum()
            mc.Add(h)

        data_fit = data.Clone(f"{data.GetName()}_fit")
        if x_range:
            data_fit = get_ranged_histogram(data_fit, x_range[0], x_range[1])

        self.data_integral = data_fit.GetSum()

        self.fitter = ROOT.TFractionFitter(data_fit, mc)
        for i, s in enumerate(self.samples):
            self.fitter.Constrain(i, 0., 1.)
            if "fixed" in self.chanenl.likelihood_fit.keys():
                if s.fitName in self.chanenl.likelihood_fit["fixed"]:
                    self.fixed += [s.fitName]
                    logger.info(f"Setting {s.name} to constant in the fit")
                    ratio = self.mc_integral[s] / self.data_integral
                    if ratio > 0:
                        self.fitter.Constrain(i, ratio * 0.9999, ratio * 1.0001)
                    else:
                        self.fitter.Constrain(i, ratio * 1.0001, ratio * 0.9999)

    def save_results(self):
        if not os.path.isdir(self.output):
            os.makedirs(self.output)
        json_dump = json.dumps(self.result)
        f = open(f"{os.path.join(self.output, self.name)}.json", "w")
        f.write(json_dump)
        f.close()

    def fit(self):
        status = self.fitter.Fit()
        status.Print()
        for i, s in enumerate(self.samples):
            frac = ROOT.Double()
            frac_err = ROOT.Double()
            self.fitter.GetResult(i, frac, frac_err)
            scale = frac * self.data_integral / self.mc_integral[s]
            scale_err = scale * frac_err / frac
            self.result[s.fitName] = (scale, scale_err)
