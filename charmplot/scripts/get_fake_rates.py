#!/usr/bin/env python
from array import array
from charmplot.common import utils
from charmplot.control import channel
from charmplot.control import sample
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


def main(options, args):
    # input file
    f = ROOT.TFile(options.input, "READ")

    # samples
    samples = options.samples.split(",")

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

        # x range and channel name
        x_range = [int(x) for x in c.split(":")[1:]]
        c = c.split(":")[0]
        print(c, x_range)

        # fake rate map
        fake_rate_map = {}
        eta_range = {}

        for s in samples:

            lep_pt_eta_tight = f.Get(f"{s}_Tight_{c}_lep_pt_eta")
            lep_pt_eta_loose = f.Get(f"{s}_AntiTight_{c}_lep_pt_eta")
            lep_pt_eta_tight.RebinX(2)
            lep_pt_eta_loose.RebinX(2)
            set_range(lep_pt_eta_tight, x_range)
            set_range(lep_pt_eta_loose, x_range)

            fake_rate = lep_pt_eta_tight.Clone(f"Fake_Rate_{s}_{c}")
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
                projX = fake_rate.ProjectionX(f"Fake_Rate_{s}_{c}_{y}", y, y)
                h_F, h_f = make_fake_rate_histograms(projX, x_range)
                fake_rate_map[y][s] = h_f

                # write to file
                h_F.Write(h_F.GetName())
                h_f.Write(h_f.GetName())

                # temporarily save in memory
                histograms[y] = h_f

            # assemble 2D histogram
            h_f_2D = make_fake_rate_histogram_2D(fake_rate, histograms)
            h_f_2D.Write()

        # plot fake rates
        for y in fake_rate_map:

            # mc_map
            mc_map = {}

            # make stack
            hs = ROOT.THStack()
            for i, s in enumerate(samples):
                samp = sample.Sample(s, None, **{'add': [s], 'subtract': []})
                fake_rate_map[y][s].SetMarkerColor(i + 1)
                fake_rate_map[y][s].SetLineColor(i + 1)
                mc_map[samp] = fake_rate_map[y][s]
                hs.Add(fake_rate_map[y][s])

            # channel object for plotting
            eta = eta_range[y]
            chan = channel.Channel(f"{s}_{c}_{y}", [c, f"|#eta| = [{eta[0]}, {eta[1]}]"], "2018", [], [])

            # canvas
            canv = utils.make_canvas(hs.GetStack().Last(), lep_pt, chan, x=800, y=800, y_split=0, events="f")
            canv.proxy_up.GetXaxis().SetRangeUser(0, x_range[1] + 20)
            canv.proxy_up.SetMinimum(y_range[0])
            canv.proxy_up.SetMaximum(y_range[1])
            canv.make_legend(None, None, mc_map, mc_map.keys(), draw_option="pe")
            hs.Draw("same nostack")
            canv.print(os.path.join(plots_folder, f"{options.output}_{c}_{y}.pdf"))

    # close out file
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
                      default="FakeRates")
    parser.add_option('-c', '--channels',
                      action="store", dest="channels",
                      help="comma separated list of channels",
                      default="2018_el_QCD_Dplus:30:90,2018_mu_QCD_Dplus:30:70")
    parser.add_option('-s', '--samples',
                      action="store", dest="samples",
                      help="comma separated list of samples",
                      default="Multijet")
    parser.add_option('-y', '--yrange',
                      action="store", dest="yrange",
                      help="y-axis range",
                      default="0:1.2")

    # parse input arguments
    options, args = parser.parse_args()

    # run
    main(options, args)
