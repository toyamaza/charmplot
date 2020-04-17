#!/usr/bin/env python
import ROOT
import os

from copy import deepcopy
from hfroofit.utils import fitUtils

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.SetAtlasStyle()


def get_channels(options):
    channels = [c.split(":")[0] for c in options.channels.split(",")]
    channel_short = dict()
    for channel in channels:
        for c in options.channels.split(","):
            if channel == c.split(":")[0]:
                if len(c.split(":")) > 1:
                    channel_short.update({channel: c.split(":")[1]})
                else:
                    channel_short.update({channel: channel})
    return channels, channel_short


def get_histograms(f, channels, channel_short, samples, var):
    histograms = dict()
    data = dict()
    for channel in channels:
        histograms.update({channel: {}})
        for sample in samples:
            h_temp = f.Get(f"{sample}_{channel}_{var}")
            h = h_temp.Clone(f"{sample}_{channel_short[channel]}_{var}")
            histograms[channel].update({sample: deepcopy(h)})
        h_temp = f.Get(f"Data_{channel}_{var}")
        h_data = h_temp.Clone(f"Data_{channel_short[channel]}_{var}")
        data.update({channel: deepcopy(h_data)})
    f.Close()
    return histograms, data


def main(options):
    # input file
    f = ROOT.TFile(options.input, "READ")

    # samples
    samples = options.samples.split(",")

    # channel
    var = options.var

    # list of all parameters
    floatPars = []

    # channels and short names
    # need short names in some cases because RooFit is not too happy about long names
    channels, channel_short = get_channels(options)

    # read input mc
    # need to deepcopy histograms so they don't get deleted when file is closed
    histograms, data = get_histograms(f, channels, channel_short, samples, var)

    # RooStats measurement
    meas = ROOT.RooStats.HistFactory.Measurement("meas", "meas")
    meas.SetOutputFilePrefix(options.name)
    meas.SetLumi(1.0)
    meas.SetLumiRelErr(0.02)
    meas.AddConstantParam("Lumi")
    meas.SetExportOnly(False)

    # Define HF samples
    Samples = dict()
    for channel in channels:
        Samples.update({channel: {}})
        for sample in samples:
            s = ROOT.RooStats.HistFactory.Sample(sample)
            print(histograms[channel][sample])
            s.SetHisto(histograms[channel][sample])
            if "Multijet" in sample:
                s.SetNormalizeByTheory(False)
                s.AddNormFactor(f"mu_QCD_{channel_short[channel]}", 1, -100, 100)
                floatPars += [f"mu_QCD_{channel_short[channel]}"]
            else:
                s.SetNormalizeByTheory(True)
            Samples[channel].update({sample: s})

    # normalisation factors
    constraints = options.constraints.split(",")
    for constraint in constraints:
        constraint = constraint.split(":")
        sample = constraint[0]
        args = [float(x) for x in constraint[1:]]
        floatPars += [f"mu_{sample}"]
        for channel in channels:
            Samples[channel][sample].AddNormFactor(f"mu_{sample}", *args)
        print(f"Set sample {sample} to {args}")

    # Define HF channels
    Channels = dict()
    for channel in channels:
        c = ROOT.RooStats.HistFactory.Channel(channel_short[channel])
        c.SetStatErrorConfig(0.05, "Poisson")
        c.SetData(data[channel])
        for sample in samples:
            c.AddSample(Samples[channel][sample])
        Channels.update({channel: c})
        meas.AddChannel(c)

    # HistoToWorkspaceFactoryFast object
    factory = ROOT.RooStats.HistFactory.HistoToWorkspaceFactoryFast(meas)

    # single channel workspaces
    results = {}
    for channel in channels:
        Channels[channel].Print()
        w = factory.MakeSingleChannelModel(meas, Channels[channel])
        res = fitUtils.run_fit(w, "obsData")
        fitUtils.results_to_json(res, floatPars, options.input, channel)
        results.update({channel: res})

    # Combined workspace
    if len(channels) > 1:
        w = factory.MakeCombinedModel(meas)
        res = fitUtils.run_fit(w, "obsData")
        fitUtils.results_to_json(res, floatPars, options.input, "combined")

    # Print results
    for channel in channels:
        print(f"{'='*10} {channel} {'='*10}")
        results[channel].Print()
    if len(channels) > 1:
        print(f"{'='*10} combined {'='*10}")
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
    parser.add_option('-f', '--fit-constraints',
                      action="store", dest="constraints",
                      help="fit constraints",
                      default="Top:1:0:10,Wjets_light_Rest:1:0:10")
    parser.add_option('-s', '--samples',
                      action="store", dest="samples",
                      help="list of samples to fit",
                      default="Wjets_light_Matched,Wjets_light_Rest,Wjets_tau,Zjets,Top,Diboson,Multijet")
    parser.add_option('-c', '--channels',
                      action="store", dest="channels",
                      help="channels to fit",
                      default="2018_el_wplusd_PT_Fit_OS:OS,2018_el_wplusd_PT_Fit_SS:SS")
    parser.add_option('-v', '--var',
                      action="store", dest="var",
                      help="variable to fit",
                      default="lep_pt")
    parser.add_option('-n', '--name',
                      action="store", dest="name",
                      help="name of the fit",
                      default="WPlusD")

    # parse input arguments
    options, args = parser.parse_args()

    # launch program
    main(options)
