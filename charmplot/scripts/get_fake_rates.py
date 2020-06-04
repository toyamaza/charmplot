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

# x range
x_range_el = [28, 100]
x_range_mu = [28, 60]


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

    # channels
    channels = options.channels.split(",")

    # out plots
    plots_folder = os.path.join(os.path.dirname(options.input), "FakeRates")
    if not os.path.isdir(plots_folder):
        os.makedirs(plots_folder)

    # out file
    out = ROOT.TFile(os.path.join(plots_folder, "FakeRates.root"), "RECREATE")

    for c in channels:

        if "_el_" in c:
            x_range = x_range_el
        else:
            x_range = x_range_mu

        lep_pt_eta_tight = f.Get(f"Multijet_{c}_lep_pt_eta")
        lep_pt_eta_loose = f.Get(f"Multijet_AntiTight_{c}_lep_pt_eta")
        # lep_pt_eta_tight.RebinX(2)
        # lep_pt_eta_loose.RebinX(2)
        set_range(lep_pt_eta_tight, x_range)
        set_range(lep_pt_eta_loose, x_range)

        fake_rate = lep_pt_eta_tight.Clone(f"Fake_Rate_{c}")
        fake_rate.Divide(lep_pt_eta_loose)
        out.cd()
        fake_rate.Write()

        for y in range(1, fake_rate.GetNbinsY() + 1):
            # eta edge
            eta = [fake_rate.GetYaxis().GetBinLowEdge(y), fake_rate.GetYaxis().GetBinLowEdge(y) + fake_rate.GetYaxis().GetBinWidth(y)]

            # channel object for plotting
            chan = channel.Channel(f"{c}_{y}", [c, f"|#eta| = [{eta[0]}, {eta[1]}]"], "2018", [], [])
            projX = fake_rate.ProjectionX(f"Fake_Rate_{c}_{y}", y, y)
            h_F, h_f = make_fake_rate_histograms(projX, x_range)

            h_F.SetLineColor(ROOT.kRed)
            h_F.SetMarkerColor(ROOT.kRed)

            # fake rate
            canv = utils.make_canvas(h_F, lep_pt, chan, x=800, y=800, y_split=0, events="f")
            canv.proxy_up.GetXaxis().SetRangeUser(0, x_range[1] + 20)
            canv.proxy_up.SetMaximum(1.2)
            canv.proxy_up.SetMinimum(0)
            # h_F.Draw("pe same")
            h_f.Draw("same")
            canv.print(os.path.join(plots_folder, f"FakeRate_{c}_{y}.pdf"))

            h_F.Write(h_F.GetName())
            h_f.Write(h_f.GetName())

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
    parser.add_option('-c', '--channels',
                      action="store", dest="channels",
                      help="comma separated list of channels",
                      default="2018_el_QCD,2018_mu_QCD")

    # parse input arguments
    options, args = parser.parse_args()

    # run
    main(options, args)
