#!/usr/bin/env python
from charmplot.common import utils
from charmplot.control import globalConfig
import os
import ROOT

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasLabels.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasUtils.C"))
ROOT.SetAtlasStyle()


def main(options, conf):

    # fitted variable
    var = conf.get_var(options.var)
    print(var.name)

    # channels
    channel = conf.get_channel(f"{options.channel}")
    print(channel.name)
    print(channel.samples)

    # label
    channel.label += ["W cuts: MET > 30, mT > 60 GeV"]

    # histos
    h_SR = {}
    h_Loose = {}

    # proxy histogram
    proxy = None

    # file
    f = ROOT.TFile(options.input, "READ")
    for s_name in channel.samples:
        s = conf.get_sample(s_name)
        h = f.Get(f"{s.shortName}_{channel.name}_{var.name}")
        if proxy is None:
            proxy = h.Clone("proxy")
        if "Loose" in s.shortName:
            h_Loose[s.shortName.replace("_Loose", "")] = h
        else:
            h_SR[s.shortName] = h

    for i in range(1, proxy.GetNbinsX() + 1):
        proxy.SetBinContent(i, 0)
        proxy.SetBinError(i, 0)

    # divide
    for s in h_SR:
        h_SR[s].Divide(h_Loose[s])
        for i in range(1, proxy.GetNbinsX() + 1):
            F = h_SR[s].GetBinContent(i)
            err = h_SR[s].GetBinError(i)
            y_ratio = F / (1 + F)
            y_ratio_err = err / (1 + F)**2
            h_SR[s].SetBinContent(i, y_ratio)
            h_SR[s].SetBinError(i, y_ratio_err)

    # mc_map
    mc_map = {}
    samples = []
    for s_name in channel.samples:
        s = conf.get_sample(s_name)
        if "Loose" not in s.shortName:
            samples += [s]
            mc_map[s] = h_SR[s.shortName]
            mc_map[s].SetLineColor(s.lineColor)
            mc_map[s].SetMarkerColor(s.lineColor)
            mc_map[s].SetMarkerStyle(24)

    # make plots
    canv = utils.make_canvas_mc_ratio(proxy, var, channel, "Ratio", x=800, y=800, events="Fraction Passing W Cuts", y_split=0.3, ratio_range=[0.89, 1.11])
    canv.pad1.cd()
    # make legend
    canv.make_legend(mc_map, samples, print_yields=False)

    # set maximum after creating legend
    canv.set_maximum([mc_map[s] for s in samples], var, mc_map[samples[0]])

    proxy.GetYaxis().SetRangeUser(0, 1.0)
    for s in mc_map:
        mc_map[s].Draw("p e same")

    # ratio plot
    h_tmp = []
    canv.pad2.cd()
    canv.proxy_dn.GetXaxis().SetNoExponent()
    # canv.proxy_dn.GetXaxis().SetMoreLogLabels()
    denum = mc_map[samples[0]]
    for i in range(0, len(mc_map)):
        h_eff_ratio = mc_map[samples[i]].Clone(f"{mc_map[samples[i]].GetName()}_ratio")
        h_eff_ratio.Divide(denum)
        h_tmp += [h_eff_ratio]
        h_eff_ratio.Draw("p e same")

    canv.print(f"{options.analysis_config.replace('.yaml', '')}/fraction_pass_w.pdf")
    canv.pad1.SetLogx()
    canv.pad2.SetLogx()
    canv.print(f"{options.analysis_config.replace('.yaml', '')}/fraction_pass_w_LOGX.pdf")


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-a', '--analysis-config',
                      action="store", dest="analysis_config",
                      help="analysis config file")
    parser.add_option('-i', '--input',
                      action="store", dest="input",
                      default="", help="Input file.")
    parser.add_option('-c', '--channel',
                      action="store", dest="channel",
                      help="channel to use",
                      default="OS-SS_Dplus")
    parser.add_option('-v', '--var',
                      action="store", dest="var",
                      help="fitted variable",
                      default="Dmeson_pt_fit")

    # parse input arguments
    options, args = parser.parse_args()

    # analysis configs
    config = options.analysis_config

    # config object
    conf = globalConfig.GlobalConfig(config, "")

    # do the plotting
    main(options, conf)
