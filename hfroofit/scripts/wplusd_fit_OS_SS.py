#!/usr/bin/env python
import ROOT
import os

from hfroofit.utils import fitUtils

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.SetAtlasStyle()

# fit variable
var = "lep_pt"

# channels
channel = "2018_el_wplusd_PT_Fit"

# list of samples
samples = [
    "Wjets_light_Matched",
    "Wjets_light_Rest",
    "Wjets_tau",
    "Zjets",
    "Top",
    "Diboson",
    "Multijet",
]


def main(options):
    # input file
    f = ROOT.TFile(options.input, "READ")

    # read input mc
    histograms_OS = {s: f.Get(f"{s}_{channel}_OS_{var}") for s in samples}
    histograms_SS = {s: f.Get(f"{s}_{channel}_SS_{var}") for s in samples}
    for s in samples:
        histograms_OS[s].SetName(f"{s}_OS")
        histograms_SS[s].SetName(f"{s}_SS")

    # read input data
    data_OS = f.Get(f"Data_{channel}_OS_{var}")
    data_SS = f.Get(f"Data_{channel}_SS_{var}")
    data_OS.SetName("Data_OS")
    data_SS.SetName("Data_SS")

    # RooStats measurement
    meas = ROOT.RooStats.HistFactory.Measurement("meas", "meas")
    meas.SetOutputFilePrefix("WPlusD")
    meas.SetLumi(1.0)
    meas.SetLumiRelErr(0.02)
    meas.AddConstantParam("Lumi")
    meas.SetExportOnly(False)
    meas.SetPOI("mu_wjets")

    # Define HF samples
    samples_OS = {}
    samples_SS = {}
    for s in samples:
        samples_OS[s] = ROOT.RooStats.HistFactory.Sample(s)
        samples_OS[s].SetNormalizeByTheory(True)
        samples_OS[s].SetHisto(histograms_OS[s])
        samples_SS[s] = ROOT.RooStats.HistFactory.Sample(s)
        samples_SS[s].SetNormalizeByTheory(True)
        samples_SS[s].SetHisto(histograms_SS[s])

    # normalisation factors
    samples_OS["Multijet"].AddNormFactor("mu_QCD_OS", 1, -1000, 1000)
    samples_OS["Wjets_light_Rest"].AddNormFactor("mu_wjets", 1, 0, 10)
    samples_OS["Top"].AddNormFactor("mu_top", 1, 0, 10)
    samples_SS["Multijet"].AddNormFactor("mu_QCD_SS", 1, -1000, 1000)
    samples_SS["Wjets_light_Rest"].AddNormFactor("mu_wjets", 1, 0, 10)
    samples_SS["Top"].AddNormFactor("mu_top", 1, 0, 10)

    # Define HF channels
    chan_OS = ROOT.RooStats.HistFactory.Channel("OS")
    chan_OS.SetStatErrorConfig(0.05, "Poisson")
    chan_OS.SetData(data_OS)
    for s in samples:
        chan_OS.AddSample(samples_OS[s])
    meas.AddChannel(chan_OS)

    chan_SS = ROOT.RooStats.HistFactory.Channel("SS")
    chan_SS.SetStatErrorConfig(0.05, "Poisson")
    chan_SS.SetData(data_SS)
    for s in samples:
        chan_SS.AddSample(samples_SS[s])
    meas.AddChannel(chan_SS)

    # HistoToWorkspaceFactoryFast object
    factory = ROOT.RooStats.HistFactory.HistoToWorkspaceFactoryFast(meas)

    # OS-only workspace
    w_OS = factory.MakeSingleChannelModel(meas, chan_OS)
    res_OS = fitUtils.run_fit(w_OS, "obsData")

    # SS-only workspace
    w_SS = factory.MakeSingleChannelModel(meas, chan_SS)
    res_SS = fitUtils.run_fit(w_SS, "obsData")

    # # Combined workspace
    w = factory.MakeCombinedModel(meas)
    res = fitUtils.run_fit(w, "obsData")

    # print results
    res_OS.Print()
    res_SS.Print()
    res.Print()


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-i', '--input',
                      action="store", dest="input",
                      help="input file with histograms")

    # parse input arguments
    options, args = parser.parse_args()

    # launch program
    main(options)
