#!/usr/bin/env python
import logging
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


def make_ratio_range(num, den, lims, name, sf):
    hist_tmp = num.Clone(name + "_tmp")
    den_Clone = den.Clone(name + "_tmp2")
    den_Clone.Scale(sf)
    hist_tmp.Divide(den_Clone)
    upper_bin = num.GetXaxis().FindBin(lims[1])
    lower_bin = num.GetXaxis().FindBin(lims[0])
    upper_bin_denom = num.GetNbinsX() - num.GetXaxis().FindBin(lims[1]) + 1
    lower_bin_denom = num.GetXaxis().FindBin(lims[0])
    nbins = upper_bin - lower_bin
    hist_out = ROOT.TH1F(name + "_final", name + "_final", nbins, lims[0], lims[1])
    for i in range(hist_tmp.GetNbinsX() - 1):
        lb_num = 0
        ub_num = 0
        lb_den = 0
        ub_den = 0
        if i + 1 < lower_bin:
            lb_num += num.GetBinContent(i + 1)
            lb_den += den_Clone.GetBinContent(i + 1)
        elif i + 1 > upper_bin:
            ub_num += num.GetBinContent(i + 1)
            ub_den += den_Clone.GetBinContent(i + 1)
        else:
            hist_out.SetBinContent(i + 2 - lower_bin, hist_tmp.GetBinContent(i + 1))
    if lb_den == 0:
        lb_den = 1
    if ub_den == 0:
        ub_den = 1
    hist_out.SetBinContent(0, (lb_num / lb_den) / lower_bin_denom)
    hist_out.SetBinContent(nbins + 1, (ub_num / ub_den) / upper_bin_denom)
    print('test')
    print(upper_bin_denom)
    return hist_out


def main(options, args):

    # input file
    f = ROOT.TFile(options.input, "READ")

    # prompt channels subtracted from Data
    prompt_channels = ["Wjets", "Zjets", "Top", "Diboson"]
    # leptons & btags
    leptons = ["mu", "el"]
    btags = ["0tag", "1tag"]
    # years and type of D
    # years = ["2016","2017", "2018"]
    # non closure calculated for 16- 18 together and 15 separately
#    years = [""]
    years = ["_2015", ""]
    # ratio of separated years to summed years to scale MJ
    lumi_ratios = [1.0, 1.0]
    dspecies = options.dspecies
    var = options.var_rw

    # out plots
    plots_folder = os.path.join(os.path.dirname(options.input), f"{options.output}")
    if not os.path.isdir(plots_folder):
        os.makedirs(plots_folder)

    # out file
    out = ROOT.TFile(os.path.join(plots_folder, f"{options.output}.root"), "RECREATE")

    data_minus_prompt_hists = []
    multijet_hists = []

    for lep in leptons:
        for t in btags:
            for y, sf in zip(years, lumi_ratios):
                # MJ hists, using 16 -18 for all years
                hist_mj = f.Get(f"Multijet_MatrixMethod{y}_{lep}_QCD_{t}_{dspecies}_{var}").Clone(f"Multijet_MatrixMethod_{lep}_QCD_{t}_{dspecies}_{var}")
                hist_mj.Sumw2()
                hist_data_minus_prompt = f.Get(f"Data{y}_{lep}_QCD_{t}_{dspecies}_{var}").Clone(f"DataMinusPrompt_{lep}_QCD_{t}_{dspecies}_{var}")
                hist_data_minus_prompt.Sumw2()
                for c in prompt_channels:
                    hist_tmp3 = f.Get(f"{c}{y}_{lep}_QCD_{t}_{dspecies}_{var}")
                    hist_data_minus_prompt.Add(hist_tmp3, -1)
                hist_data_minus_prompt.Write()
                data_minus_prompt_hists += [hist_data_minus_prompt]
                # MET RW histograms
                ratio = hist_data_minus_prompt.Clone(f"RW_ratio_{lep}_QCD_{t}_{dspecies}_{var}")
                ratio.Sumw2()
                ratio.Divide(hist_mj)
                # ratio.Scale(1./sf)
                ratio.Write()
                make_ratio_range(hist_data_minus_prompt, hist_mj, options.var_rw_range, f"RW_ratio{y}_{lep}_QCD_{t}_{dspecies}_{var}", sf).Write()
                multijet_hists += [hist_mj]
                hist_mj.Write()
    out.Close()


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-i', '--input',
                      action="store", dest="input",
                      help="input root file")
    parser.add_option('-o', '--output',
                      action="store", dest="output",
                      help="output name",
                      default="MMReweight")
    parser.add_option('-d', '--dspecies',
                      action="store", dest="dspecies",
                      help="which d meson",
                      default="Dplus")
    parser.add_option('-v', '--variable-reweight',
                      action="store", dest="var_rw",
                      help="variable to reweight wrt (MET default)",
                      default="met_met")
    parser.add_option('-r', '--variable-reweight-range',
                      action="store", dest="var_rw_range",
                      help="python list range ofvariable to reweight wrt (MET default)",
                      default=[0, 60])

    # parse input arguments
    options, args = parser.parse_args()

    # config
    # from charmplot.control import globalConfig
    # conf = globalConfig.GlobalConfig(options.analysis_config, options.analysis_config)

    # run
    main(options, args)
