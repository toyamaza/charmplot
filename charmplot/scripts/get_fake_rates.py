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


def make_fake_rate_histogram_2D(h, histograms):
    xbins = []
    for i in range(1, histograms[1].GetNbinsX() + 1):
        xbins += [histograms[1].GetBinLowEdge(i)]
    xbins += [histograms[1].GetBinLowEdge(len(xbins) + 1)]

    ybins = []
    for i in range(1, h.GetNbinsY() + 1):
        ybins += [h.GetYaxis().GetBinLowEdge(i)]
    ybins += [h.GetYaxis().GetBinLowEdge(len(ybins) + 1)]

    hout = ROOT.TH2D(f"{h.GetName()}_final", "", len(xbins) - 1, array('d', xbins), len(ybins) - 1, array('d', ybins))
    for i in range(1, h.GetNbinsY() + 1):
        for j in range(0, histograms[i].GetNbinsX() + 2):
            hout.SetBinContent(j, i, histograms[i].GetBinContent(j))
            hout.SetBinError(j, i, histograms[i].GetBinError(j))

    return hout


def make_fake_rate_histograms(h, x_range):
    # even bins
    if h.GetBinWidth(1) == h.GetBinWidth(h.GetNbinsX()):
        bin_width = h.GetBinWidth(1)
        nbins = (x_range[1] - x_range[0]) / bin_width
        h_F = ROOT.TH1F(f"{h.GetName()}_F", "", int(nbins), x_range[0], x_range[1])
        h_f = ROOT.TH1F(f"{h.GetName()}_f", "", int(nbins), x_range[0], x_range[1])
    
    # uneven bins
    else:
        xbins = []
        for i in range(1, h.GetNbinsX() + 1):
            xbins += [h.GetBinLowEdge(i)]
        xbins += [h.GetBinLowEdge(len(xbins) + 1)]
        nbins = len(xbins) - 1
        h_F = ROOT.TH1F(f"{h.GetName()}_F", "", int(nbins), array('d', xbins))
        h_f = ROOT.TH1F(f"{h.GetName()}_f", "", int(nbins), array('d', xbins))

    for i in range(1, h.GetNbinsX() + 1):
        center = h.GetBinCenter(i)
        bin_f = h_f.FindBin(center)
        if bin_f > 0 and bin_f <= h_f.GetNbinsX():
            val = h.GetBinContent(i)
            err = h.GetBinError(i)
            h_F.SetBinContent(bin_f, val)
            h_F.SetBinError(bin_f, err)
            f_val = (1 + val)
            #if f_val < 0.0:
            #    print("look here")
            #    f_val = 0.0 
            h_f.SetBinContent(bin_f, val / f_val)
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


def make_error_band(h_nom, hists):
    if len(hists) < 1:
        print("ERROR: Error band maker needs non empty list")
        return
    
    # Clone one hist to use as output template
    h_out_up = h_nom.Clone(h_nom.GetName()+"_GEN_sys_UP")
    h_out_dn = h_nom.Clone(h_nom.GetName()+"_GEN_sys_DN")
    gr_out_to_plot = ROOT.TGraphAsymmErrors()
    # find max in each bin of all histograms in the list
    first_hist = True
    for h in hists:
        for y_ind in range(1, h_nom.GetNbinsY()+1):
            for x_ind in range(1, h_nom.GetNbinsX()+1):
                base_val = h_nom.GetBinContent(x_ind, y_ind)
                if h.GetBinContent(x_ind, y_ind) < 0.0:
                    h.SetBinContent(x_ind, y_ind, 0.0)
                if h.GetBinContent(x_ind, y_ind) - base_val > h_out_up.GetBinContent(x_ind, y_ind) - base_val :
                    if h.GetBinContent(x_ind, y_ind) > 1:
                        h_out_up.SetBinContent(x_ind, y_ind, 1.0)
                    else:
                        h_out_up.SetBinContent(x_ind, y_ind, h.GetBinContent(x_ind, y_ind))
                if base_val - h.GetBinContent(x_ind, y_ind) > base_val - h_out_dn.GetBinContent(x_ind, y_ind) :
                    h_out_dn.SetBinContent(x_ind, y_ind, h.GetBinContent(x_ind, y_ind))
                
                # Dont want plotting hist for 2D case, only variations
                if h_nom.GetNbinsY() > 1: 
                    continue
                if first_hist:
                    gr_out_to_plot.SetPoint(gr_out_to_plot.GetN(), h_nom.GetBinCenter(x_ind), h_nom.GetBinContent(x_ind))
                x_err = h_nom.GetBinWidth(x_ind)/2.
                gr_out_to_plot.SetPointError(x_ind-1, x_err, x_err, h_nom.GetBinContent(x_ind) - h_out_dn.GetBinContent(x_ind), h_out_up.GetBinContent(x_ind) - h_nom.GetBinContent(x_ind))

        first_hist = False
    print("test ", f"{h_out_dn.GetNbinsX()} : {gr_out_to_plot.GetN()} ")
    return h_out_dn, h_out_up, gr_out_to_plot


def process_rebin(h_in, rebin):
    try:
        rebin_int = int(rebin)
        return h_in.RebinX(rebin_int)
    except ValueError:
        NbinsY = h_in.GetNbinsY()
        xbins = [int(x) for x in rebin.split(";")]
        NbinsX = len(xbins)

        ybins = []
        for i in range(1, h_in.GetNbinsY() + 1):
            ybins += [h_in.GetYaxis().GetBinLowEdge(i)]
        ybins += [h_in.GetYaxis().GetBinLowEdge(len(ybins) + 1)]

        hout = ROOT.TH2D(f"{h_in.GetName()}", "", len(xbins) - 1, array('d', xbins), len(ybins) - 1, array('d', ybins))
        for j in range(1, hout.GetNbinsY() + 1):
            for i in range(0,len(xbins)):
                bin_err = ROOT.Double()
                x_bin_start = xbins[i]+1
                if i < len(xbins)-1:
                    x_bin_end = xbins[i+1]
                else:
                    x_bin_end = -1
                bin_val = h_in.IntegralAndError(x_bin_start, x_bin_end, j, j, bin_err)
                hout.SetBinContent(i+1, j, bin_val)
                hout.SetBinError(i+1, j, bin_err)
            # rightmost overflow bin
            #bin_err = ROOT.Double()
            #bin_val = h_in.IntegralAndError(xbins[-1], -1, j, j, bin_err)
            #hout.SetBinContent(i+1, j, bin_val)
            #hout.SetBinError(i+1, j, bin_err)
        print("here", hout.GetBinContent(3,1))
        return hout


def main(conf, options, args):
    # input file
    f = ROOT.TFile(options.input, "READ")

    # samples
    samples = options.samples.split(",")
    samples_base = samples[0]
    
    #systematics
    systematics = options.systematics.split(",")
    if systematics[0] != "":
        systematics.insert(0,"")

    # channels
    channels = options.channels.split(",")

    # out plots
    plots_folder = os.path.join(os.path.dirname(options.input), f"{options.output}")
    if not os.path.isdir(plots_folder):
        os.makedirs(plots_folder)

    # out file
    out = ROOT.TFile(os.path.join(plots_folder, f"{options.output}.root"), "RECREATE")

    # y range
    y_range = [float(x) for x in options.yrange.split(":")]

    for c in channels:

        # channel name
        channel_name = c.split(":")[0]

        # x range and channel name
        x_range = [int(x) for x in c.split(":")[1:3]]

        # fake rate map
        fake_rate_map = {}
        eta_range = {}

        for s in samples:
            h_sys = []
            for sys in systematics:
                if s != samples_base and sys != "":
                    continue
                if sys != "":
                    sys = "_"+sys
                v = "lep_pt_eta" if "el_" not in channel_name else "lep_pt_calo_eta"
                # v = "lep_pt_met" if "el_" not in channel_name else "lep_pt_met"
                print(f"{s}_Tight_{channel_name}_{v}{sys}")
                lep_pt_eta_tight = f.Get(f"{s}_Tight_{channel_name}_{v}{sys}")
                lep_pt_eta_loose = f.Get(f"{s}_AntiTight_{channel_name}_{v}{sys}")
                rebin = c.split(":")[3]
                #rebin = int(c.split(":")[3])
                lep_pt_eta_tight = process_rebin(lep_pt_eta_tight, rebin)
                lep_pt_eta_loose = process_rebin(lep_pt_eta_loose, rebin)
                #lep_pt_eta_tight.RebinX(rebin)
                #lep_pt_eta_loose.RebinX(rebin)
                set_range(lep_pt_eta_tight, x_range)
                set_range(lep_pt_eta_loose, x_range)
                fake_rate = lep_pt_eta_tight.Clone(f"Fake_Rate_{s}{sys}_{channel_name}")
                fake_rate.Divide(lep_pt_eta_loose)
                out.cd()
                fake_rate.Write()
                # store fake factor histograms for later
                histograms = {}
                for y in range(1, fake_rate.GetNbinsY() + 1):
                    # insert y slice in map
                    if y not in fake_rate_map.keys():
                        fake_rate_map[y] = {}
                    # eta edge
                    eta = [fake_rate.GetYaxis().GetBinLowEdge(y), fake_rate.GetYaxis().GetBinLowEdge(y) + fake_rate.GetYaxis().GetBinWidth(y)]
                    eta_range[y] = eta
                    # get fake rate
                    projX = fake_rate.ProjectionX(f"Fake_Rate_{s}{sys}_{channel_name}_{y}", y, y)
                    h_F, h_f = make_fake_rate_histograms(projX, x_range)
                    fake_rate_map[y][s+sys] = h_f
                    # write to file
                    h_F.Write(h_F.GetName())
                    h_f.Write(h_f.GetName())
                    # temporarily save in memory
                    histograms[y] = h_f
                # assemble 2D histogram
                h_f_2D = make_fake_rate_histogram_2D(fake_rate, histograms)
                if sys != "":
                    h_sys += [h_f_2D]
                else:
                    h_f_2D.Write()

                if s == samples_base and sys == "":
                    h_nom_tmp = h_f_2D

            if s == samples_base and len(systematics) > 1:
                h_sys_band_dn, h_sys_band_up, gr_sys_band_to_plot = make_error_band(h_nom_tmp, h_sys)
                out.cd()
                h_sys_band_dn.Write()
                h_sys_band_up.Write()

        # plot fake rates
        for y in fake_rate_map:

            # mc_map
            mc_map = {}

            # make stack
            hs = ROOT.THStack()
            # make systematic list
            h_sys = []
            for s in samples:
                for sys in systematics:
                    if s != samples_base or sys == "":
                        continue
                    if sys != "":
                        sys = "_"+sys
                    h_sys += [fake_rate_map[y][s+sys]]
                sample_config = conf.get_sample(s)
                fake_rate_map[y][s].SetMarkerColor(sample_config.lineColor)
                fake_rate_map[y][s].SetLineColor(sample_config.lineColor)
                if s == samples_base and len(systematics) > 1:
                    h_sys_band_dn, h_sys_band_up, gr_sys_band_to_plot = make_error_band(fake_rate_map[y][s], h_sys)
                    h_sys_band_dn.Write()
                    h_sys_band_up.Write()
                mc_map[sample_config] = fake_rate_map[y][s]
                hs.Add(fake_rate_map[y][s])
                print(f"y: {y}, s: {s}, value of last bin{fake_rate_map[y][s].GetBinContent(fake_rate_map[y][s].GetNbinsX()+1)}")

            # channel object for plotting
            eta = eta_range[y]
            chan = channel.Channel(f"{s}_{channel_name}_{y}", [channel_name, f"|#eta| = [{eta[0]}, {eta[1]}]"], "2016+2017+2018", [], [])

            # canvas
            var = lep_pt
            canv = utils.make_canvas(hs.GetStack().Last(), var, chan, x=800, y=800, y_split=0, events="f")
            canv.proxy_up.GetXaxis().SetRangeUser(0, x_range[1] + 20)
            canv.proxy_up.SetMinimum(y_range[0])
            canv.proxy_up.SetMaximum(y_range[1])
            canv.make_legend(None, None, mc_map, mc_map.keys(), draw_option="pe", leg_offset=-0.05)
            #h_sys_band.SetLineColor(1)
            #_sys_band.Draw("same")
            if len(systematics) > 1:
                gr_sys_band_to_plot.SetFillColor(1)
                gr_sys_band_to_plot.SetFillStyle(3354)
                gr_sys_band_to_plot.Draw("e2")
            #gr_sys_band_to_plot.Write()
            hs.Draw("same nostack")


            canv.print(os.path.join(plots_folder, f"{options.output}_{channel_name}_{y}.pdf"))
            canv.print(os.path.join(plots_folder, f"{options.output}_{channel_name}_{y}.png"))
    
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
    parser.add_option('--systematics',
                      action="store", dest="systematics",
                      help="comma separated list of systematics, to be applied to baseline",
                      default="")
    parser.add_option('-y', '--yrange',
                      action="store", dest="yrange",
                      help="y-axis range",
                      default="0:1.2")

    # parse input arguments
    options, args = parser.parse_args()

    # config
    from charmplot.control import globalConfig
    conf = globalConfig.GlobalConfig(options.analysis_config, options.analysis_config)

    # run
    main(conf, options, args)
