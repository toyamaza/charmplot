#!/usr/bin/env python
from array import array
from charmplot.common import utils
from charmplot.control import channel
from charmplot.control import variable
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

# lepton pt variable
lep_pt = variable.Variable("lep_pt", **{
    "label": "p_{T}",
    "unit": "GeV",
})

# horizontal error bars for histograms
ROOT.gStyle.SetErrorX(0.5)


def make_fake_rate_histogram_2D(h, histograms, ybins, xbins):

    for i in range(1, histograms[1].GetNbinsX() + 1):
        xbins += [histograms[1].GetBinLowEdge(i)]
    xbins += [histograms[1].GetBinLowEdge(len(xbins) + 1)]

    for i in range(1, h.GetNbinsY() + 1):
        ybins += [h.GetYaxis().GetBinLowEdge(i)]
    ybins += [h.GetYaxis().GetBinLowEdge(len(ybins) + 1)]
    hout = ROOT.TH2D(f"{h.GetName()}_final", "", len(xbins) - 1, array('d', xbins), len(ybins) - 1, array('d', ybins))
    for i in range(1, h.GetNbinsY() + 1):
        for j in range(0, histograms[i].GetNbinsX() + 2):
            hout.SetBinContent(j, i, histograms[i].GetBinContent(j))
            hout.SetBinError(j, i, histograms[i].GetBinError(j))

    return hout


def make_fake_rate_histogram_2D_TEff(eff, xbins, ybins, name, tag_index):

    graphs = []
    if tag_index < 0:
        first_eff = True
        for i in range(len(eff)):
            for j in range(len(eff[0])):
                if first_eff:
                    eff_tmp = eff[i][j].Clone(f"new_name_{i}_{j}")
                else:
                    eff_tmp += eff[i][j]
            graphs += [eff_tmp.CreateGraph()]
    else:
        for i in range(len(eff)):
            graphs += [eff[i][tag_index].CreateGraph()]

    hout = ROOT.TH2D(name + "_final", "", len(xbins) - 1, array('d', xbins), len(ybins) - 1, array('d', ybins))

    for i in range(0, len(ybins) - 1):
        # if graphs[i].GetN() < 1:
        #    print(i, ybins[i])
        for j in range(0, len(xbins) - 1):
            if graphs[i].GetN() < 1:
                hout.SetBinContent(j + 1, i + 1, 0)
            else:
                hout.SetBinContent(j + 1, i + 1, graphs[i].GetY()[j])
                hout.SetBinError(j + 1, i + 1, 0.5 * (graphs[i].GetErrorYhigh(j) + graphs[i].GetErrorYlow(j)))
        hout.SetBinContent(0, i + 1, hout.GetBinContent(1, i + 1))
        hout.SetBinContent(hout.GetNbinsX() + 1, i + 1, hout.GetBinContent(hout.GetNbinsX(), i + 1))
        hout.SetBinError(0, i + 1, hout.GetBinError(1, i + 1))
        hout.SetBinError(hout.GetNbinsX() + 1, i + 1, hout.GetBinError(hout.GetNbinsX(), i + 1))

    for j in range(0, len(xbins) - 1):
        hout.SetBinContent(j + 1, 0, hout.GetBinContent(j + 1, 1))
        hout.SetBinContent(j + 1, hout.GetNbinsY() + 1, hout.GetBinContent(j + 1, hout.GetNbinsY()))
        hout.SetBinError(j + 1, 0, hout.GetBinError(j + 1, 1))
        hout.SetBinError(j + 1, hout.GetNbinsY() + 1, hout.GetBinError(j + 1, hout.GetNbinsY()))

    hout.SetBinContent(0, 0, hout.GetBinContent(1, 1))
    hout.SetBinContent(0, hout.GetNbinsY() + 1, hout.GetBinContent(1, hout.GetNbinsY()))
    hout.SetBinError(0, 0, hout.GetBinError(1, 1))
    hout.SetBinError(0, hout.GetNbinsY() + 1, hout.GetBinError(1, hout.GetNbinsY()))
    hout.SetBinContent(hout.GetNbinsX() + 1, 0, hout.GetBinContent(hout.GetNbinsX(), 1))
    hout.SetBinContent(hout.GetNbinsX() + 1, hout.GetNbinsY() + 1, hout.GetBinContent(hout.GetNbinsX(), hout.GetNbinsY()))
    hout.SetBinError(hout.GetNbinsX() + 1, 0, hout.GetBinError(hout.GetNbinsX(), 1))
    hout.SetBinError(hout.GetNbinsX() + 1, hout.GetNbinsY() + 1, hout.GetBinError(hout.GetNbinsX(), hout.GetNbinsY()))

    return hout


def make_fake_rate_histograms(h, x_range):
    bin_width = h.GetBinWidth(1)
    nbins = (x_range[1] - x_range[0]) / bin_width
    h_F = ROOT.TH1F(f"{h.GetName()}_F", "", int(nbins), x_range[0], x_range[1])
    h_f = ROOT.TH1F(f"{h.GetName()}_f", "", int(nbins), x_range[0], x_range[1])

    for i in range(1, h.GetNbinsX() + 1):
        center = h.GetBinCenter(i)
        bin_f = h_f.FindBin(center)
        if bin_f > 0 and bin_f <= h_f.GetNbinsX():
            val = h.GetBinContent(i)
            err = h.GetBinError(i)
            h_F.SetBinContent(bin_f, val)
            h_F.SetBinError(bin_f, err)
            h_f.SetBinContent(bin_f, val / (1 + val))
            h_f.SetBinError(bin_f, err / (1 + val)**2)

    h_f.SetBinContent(0, h_f.GetBinContent(1))
    h_f.SetBinError(0, h_f.GetBinError(1))
    h_f.SetBinContent(h_f.GetNbinsX() + 1, h_f.GetBinContent(h_f.GetNbinsX()))
    h_f.SetBinError(h_f.GetNbinsX() + 1, h_f.GetBinError(h_f.GetNbinsX()))

    h_F.SetBinContent(0, h_F.GetBinContent(1))
    h_F.SetBinError(0, h_F.GetBinError(1))
    h_F.SetBinContent(h_F.GetNbinsX() + 1, h_F.GetBinContent(h_F.GetNbinsX()))
    h_F.SetBinError(h_F.GetNbinsX() + 1, h_F.GetBinError(h_F.GetNbinsX()))

    return h_F, h_f


def get_div_graphs(eff_main, eff_list):
    den = ROOT.TGraphAsymmErrors()
    den = eff_main.CreateGraph()
    out_graphs = []
    print(den.GetN())
    for num_eff in eff_list:
        num = ROOT.TGraphAsymmErrors()
        num = num_eff.CreateGraph()
        quotient = ROOT.TGraphAsymmErrors(den.GetN())
        if num.GetN() > 1:
            print("num.GetN")
            print(num.GetN())
            print("den.GetN")
            print(den.GetN())
            for i in range(min(den.GetN(),num.GetN())):
                quotient.SetPoint(i, den.GetX()[i], num.GetY()[i] / den.GetY()[i])
                quotient.SetPointError(i, den.GetErrorXlow(i), den.GetErrorXhigh(i), num.GetErrorYlow(i) / num.GetY()[i], num.GetErrorYhigh(i) / num.GetY()[i])
            quotient.SetMarkerColor(num.GetMarkerColor())
            quotient.SetLineColor(num.GetLineColor())
        out_graphs += [quotient]
    return out_graphs


def set_range(h, x_range):
    bin1 = 1
    while h.GetXaxis().GetBinCenter(bin1) < x_range[1]:
        bin1 += 1
    bin2 = h.GetNbinsX() + 1

    print(bin1, h.GetXaxis().GetBinCenter(bin1))

    for i in range(1, h.GetNbinsY() + 1):
        err = ROOT.Double()
        integral = h.IntegralAndError(bin1, bin2, i, i, err)
        err_last = ROOT.Double()
        integral_last = h.IntegralAndError(bin1 - 1, bin1 - 1, i, i, err_last)
        h.SetBinContent(bin1, i, integral + integral_last)
        h.SetBinError(bin1, i, (err**2 + err_last**2)**(0.5))
        for j in range(bin1, bin2):
            h.SetBinContent(j, i, 0)
            h.SetBinError(j, i, 0)


def main(conf, options, args):
    # input file
    f = ROOT.TFile(options.input, "READ")

    #systematic variations
    samples_sys = options.samples_sys.split(",")

    # samples
    samples = options.samples.split(",")
    # samples must include systematics
    samples += samples_sys

    # channels
    channels = options.channels.split(",")

    # d species
    dspecies = options.dspecies.split(",")

    # split tags
    tags = options.tags.split(",")

    # out plots
    plots_folder = os.path.join(os.path.dirname(options.input), f"{options.output}")
    if not os.path.isdir(plots_folder):
        os.makedirs(plots_folder)

    # out file
    out = ROOT.TFile(os.path.join(plots_folder, f"{options.output}.root"), "RECREATE")

    # y range
    y_range = [float(x) for x in options.yrange.split(":")]

    for c in channels:

        # to capture things once per channel
        once_per_channel = True

        # channel name
        channel_name = c.split(":")[0]

        # x range and channel name
        x_range = [int(x) for x in c.split(":")[1:3]]

        # fake rate map
        fake_rate_map_hists = {}
        fake_rate_map = {}
        # Tefficiencies to skip for weight reasons
        fake_rate_map_skip = {}
        eta_range = {}

        for s in samples:
            for d in dspecies:
                for t in tags:
                    v = "lep_pt_eta" if "el_" not in channel_name else "lep_pt_calo_eta"

                    rebin = int(c.split(":")[3])

                    print(f"{s}_Tight_{channel_name}_{t}_{d}_{v}")
                    lep_pt_eta_tight = f.Get(f"{s}_Tight_{channel_name}_{t}_{d}_{v}")
                    lep_pt_eta_loose = f.Get(f"{s}_AntiTight_{channel_name}_{t}_{d}_{v}")
                    lep_pt_eta_tight.RebinX(rebin)
                    lep_pt_eta_loose.RebinX(rebin)
                    set_range(lep_pt_eta_tight, x_range)
                    set_range(lep_pt_eta_loose, x_range)

                    fake_rate = lep_pt_eta_tight.Clone(f"Fake_Rate_{s}_{channel_name}_{t}_{d}")
                    fake_rate.Divide(lep_pt_eta_loose)
                    out.cd()
                    fake_rate.Write()

                    # store fake factor histograms for later
                    histograms_f = {}
                    effs = {}

                    for y in range(1, fake_rate.GetNbinsY() + 1):
                        # insert y slice in map
                        if y not in fake_rate_map.keys():
                            fake_rate_map[y] = {}
                            fake_rate_map_hists[y] = {}
                            fake_rate_map_skip[y] = {}

                        # eta edge
                        eta = [fake_rate.GetYaxis().GetBinLowEdge(y), fake_rate.GetYaxis().GetBinLowEdge(y) + fake_rate.GetYaxis().GetBinWidth(y)]
                        eta_range[y] = eta

                        # get fake rate
                        projX = fake_rate.ProjectionX(f"Fake_Rate_{s}_{channel_name}_{t}_{d}_{y}", y, y)
                        h_F, h_f = make_fake_rate_histograms(projX, x_range)
                        fake_rate_map_hists[y][f"{s}_{t}_{d}"] = h_f

                        # get fake rate via TEfficiency
                        projX_num = lep_pt_eta_tight.ProjectionX(f"Eff_num_{s}_{channel_name}_{t}_{d}_{y}", y, y)
                        projX_den = lep_pt_eta_loose.ProjectionX(f"Eff_den_{s}_{channel_name}_{t}_{d}_{y}", y, y)
                        projX_den.Add(projX_num)
                        print("before" + str(y))
                        print(projX_num.GetNbinsX())
                        print(projX_den.GetNbinsX())
                        tmp_out = projX_num.Clone(f"tmp_out_{s}_{channel_name}_{t}_{d}_{y}")
                        tmp_out.Divide(projX_den)
                        eff = ROOT.TEfficiency(projX_num, projX_den)
                        fake_rate_map[y][f"{s}_{t}_{d}"] = eff
                        # skip Tefficiency that fail construction for stats reasons
                        if eff.CheckConsistency(projX_num, projX_den):
                            fake_rate_map_skip[y][f"{s}_{t}_{d}"] = False
                        else:
                            fake_rate_map_skip[y][f"{s}_{t}_{d}"] = True
                        print("after")

                        # write to file
                        h_F.Write(h_F.GetName())
                        h_f.Write(h_f.GetName())
                        eff.Write(eff.GetName())
                        tmp_out.Write(tmp_out.GetName())

                        # temporarily save in memory
                        histograms_f[y] = h_f
                        effs[y] = eff

                    # assemble 2D histogram
                    if once_per_channel:
                        ybins_eta_chan = []
                        xbins_pt_chan = []
                        h_f_2D = make_fake_rate_histogram_2D(fake_rate, histograms_f, ybins_eta_chan, xbins_pt_chan)
                    else:
                        once_per_channel = False
                        ybins_eta = []
                        xbins_pt = []
                        h_f_2D = make_fake_rate_histogram_2D(fake_rate, histograms_f, ybins_eta, xbins_pt)

                    h_f_2D.Write()

        # plot fake rates

        eff_sums = []
        eff_sums_sys = []
        for y in fake_rate_map:

            # mc_map
            mc_map = {}

            # make stack
            fake_rate_map_comb = {}

            eff_sums_tmp = []
            eff_sums_tmp_sys = []

            for t in tags:
                first_sample = True
                first_sys_sample = True
                # eff_sum = ROOT.TEfficiency()
                for s in samples:
                    print(f"look here: {s}_{samples_sys}")
                    first_d_species = True
                    # eff_tmp = ROOT.TEfficiency()
                    for d in dspecies:
                        if fake_rate_map_skip[y][f"{s}_{t}_{d}"]:
                            continue
                        if first_d_species:
                            eff_tmp = fake_rate_map[y][f"{s}_{t}_{d}"].Clone(f"{s}_{t}_{d}_Comb")
                            print("first D species: " + f"{s}_{t}_{d}")
                            first_d_species = False
                        else:
                            eff_tmp.Add(fake_rate_map[y][f"{s}_{t}_{d}"])
                            print("second D species: " + f"{s}_{t}_{d}")
                    sample_config = conf.get_sample(f"{s}_{t}")
                    fake_rate_map_comb[f"{s}_{t}"] = eff_tmp
                    fake_rate_map_comb[f"{s}_{t}"].SetMarkerColor(sample_config.lineColor)
                    fake_rate_map_comb[f"{s}_{t}"].SetLineColor(sample_config.lineColor)
                    mc_map[sample_config] = fake_rate_map_comb[f"{s}_{t}"]

                    # skip the systematic variation in the average, make own avg
                    if s in samples_sys:
                        if first_sys_sample:
                            eff_sum_sys = eff_tmp.Clone(f"{s}_{t}_Comb")
                            first_sys_sample = False
                        else:
                            eff_sum_sys.Add(eff_tmp)
                    else:
                        # Prepping sum for output
                        if first_sample:
                            eff_sum = eff_tmp.Clone(f"{s}_{t}_Comb")
                            first_sample = False
                        else:
                            eff_sum.Add(eff_tmp)
                # for all tags, summing efficiences
                eff_sums_tmp += [eff_sum]
                eff_sums_tmp_sys += [eff_sum_sys]
            # add all tags into eta list to be exported
            eff_sums += [eff_sums_tmp]
            eff_sums_sys += [eff_sums_tmp_sys]

            sample_config = conf.get_sample("MC_Comb")
            eff_sum.SetMarkerColor(sample_config.lineColor)
            eff_sum.SetLineColor(sample_config.lineColor)
            mc_map[sample_config] = eff_sum

            # channel object for plotting
            eta = eta_range[y]
            chan = channel.Channel(f"{s}_{channel_name}_{y}", [channel_name, f"|#eta| = [{eta[0]}, {eta[1]}]"], "2016+2017+2018", [], [])

            # canvas
            var = lep_pt
            canv = utils.make_canvas(histograms_f[y], var, chan, x=800, y=800, events="f")
            canv.proxy_up.GetXaxis().SetRangeUser(0, x_range[1] + 20)
            canv.proxy_dn.GetXaxis().SetRangeUser(0, x_range[1] + 20)
            canv.proxy_up.SetMinimum(y_range[0])
            canv.proxy_up.SetMaximum(y_range[1])
            canv.make_legend(None, None, mc_map, mc_map.keys(), draw_option="pe", leg_offset=-0.15)

            eff_list = []
            for s in samples:
                for t in tags:
                    fake_rate_map_comb[f"{s}_{t}"].Draw("same")
                    eff_list += [fake_rate_map_comb[f"{s}_{t}"]]
            eff_sum.Draw("same")
            # Draw 2nd pane
            canv.pad2.cd()
            canv.set_ratio_range(.95, 1.05, override=True)
            canv.proxy_dn.GetYaxis().SetTitle("Sample/Comb")

            graph_list = get_div_graphs(eff_sum, eff_list)
            for gr in graph_list:
                gr.Draw("psame")
                # h_ratio = utils.make_ratio(fake_rate_map_hists[y][f"{s}_{t}_{d}"], fake_rate_map_hists[y][f"{samples[0]}_{dspecies[0]}"]).Clone()
                # h_ratios += [h_ratio]
                # h_ratio.Draw("same")

            canv.print(os.path.join(plots_folder, f"{options.output}_{channel_name}_{y}.pdf"))
            canv.print(os.path.join(plots_folder, f"{options.output}_{channel_name}_{y}.png"))

        first_tag = True
        for t, i in zip(tags, range(len(tags))):
            hf_2D_eff = make_fake_rate_histogram_2D_TEff(eff_sums, xbins_pt_chan, ybins_eta_chan, f"Fake_Rate_{channel_name}_{t}", i)
            hf_2D_eff.Write()
            hf_2D_eff_sys = make_fake_rate_histogram_2D_TEff(eff_sums_sys, xbins_pt_chan, ybins_eta_chan, f"Fake_Rate_{channel_name}_{t}_sys", i)
            hf_2D_eff_sys.Write()

        eff_tag_sum = make_fake_rate_histogram_2D_TEff(eff_sums, xbins_pt_chan, ybins_eta_chan, f"Fake_Rate_{channel_name}", -1)
        eff_tag_sum_sys = make_fake_rate_histogram_2D_TEff(eff_sums_sys, xbins_pt_chan, ybins_eta_chan, f"Fake_Rate_{channel_name}_sys", -1)
        eff_tag_sum.Write()
        eff_tag_sum_sys.Write()

    # close out file
    out.Close()


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
                      help="input root file")
    parser.add_option('-o', '--output',
                      action="store", dest="output",
                      help="output name",
                      default="FakeRates")
    parser.add_option('-c', '--channels',
                      action="store", dest="channels",
                      help="comma separated list of channels",
                      default="el_QCD_0tag_Dplus:30:90,el_QCD_1tag_Dplus:30:90,mu_QCD_0tag_Dplus:30:70,mu_QCD_1tag_Dplus:30:70")
    parser.add_option('-s', '--samples',
                      action="store", dest="samples",
                      help="comma separated list of samples",
                      default="Multijet")
    parser.add_option('-y', '--yrange',
                      action="store", dest="yrange",
                      help="y-axis range",
                      default="0:1.2")
    parser.add_option('-d', '--dspecies',
                      action="store", dest="dspecies",
                      help="comma separated list D species to be compared",
                      default="Dplus")
    parser.add_option('-t', '--tags',
                      action="store", dest="tags",
                      help="comma separated list tags to be compared",
                      default="0tag")
    parser.add_option('-k', '--samples_sys',
                      action="store", dest="samples_sys",
                      help="MC that will be used as systematic variation, not to be included in general sum",
                      default="MG_Wjets_emu")

    # parse input arguments
    options, args = parser.parse_args()

    # config
    from charmplot.control import globalConfig
    conf = globalConfig.GlobalConfig(options.analysis_config, options.analysis_config)

    # run
    main(conf, options, args)
