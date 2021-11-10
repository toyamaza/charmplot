#!/usr/bin/env python
import ROOT
import os
dirname = os.path.join(os.path.dirname(__file__), "../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasLabels.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasUtils.C"))
ROOT.SetAtlasStyle()


def createCanvasPads():
    c = ROOT.TCanvas("c", "canvas", 800, 800)
    # Upper histogram plot is pad1
    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.2)
    pad2.Draw()

    return c, pad1, pad2


def main(options, args):

    f = ROOT.TFile(os.path.join(options.input, "histograms.root"))

    h1 = f.Get(f"{options.h1}").Clone(f"{options.h1}")
    h2 = f.Get(f"{options.h2}").Clone(f"{options.h2}")

    h2.Scale(h1.GetSumOfWeights() / h2.GetSumOfWeights())
    h1.Scale(1 / h1.GetSumOfWeights())
    h2.Scale(1 / h2.GetSumOfWeights())

    h1.SetFillStyle(0)
    h1.SetLineColor(ROOT.kRed)
    h2.SetFillStyle(0)
    h2.SetLineColor(ROOT.kBlue)

    c, pad1, pad2 = createCanvasPads()
    pad1.cd()
    h1.Draw("hist E")
    h2.Draw("hist E same")
    ROOT.ATLASLabel(0.17, 0.90, "Internal", 1)
    ROOT.myText(0.17, 0.9 - 1 * 0.06, 1, "#sqrt{s} = 13 TeV")

    comp_name = (options.h1.split("_0tag")[0]) + "_" + (options.h2.split("_0tag")[0]).split("OS_")[-1] + "_" + (options.h1.split("_Dmeson")[0]).split("s_")[-1]
    ROOT.myText(0.17, 0.9 - 2 * 0.06, 1, comp_name.split("_OS")[0])
    ROOT.myText(0.17, 0.9 - 3 * 0.06, 1, (comp_name.split("_0tag")[0]).split("OS_")[-1])

    pad2.cd()
    hr = h2.Clone("hr")
    hr.Add(h1, -1)
    hr.Divide(h1)
    [hr.AddBinContent(i) for i in range(1, hr.GetNbinsX() + 1)]
    hr.Draw()
    hr.GetYaxis().SetRangeUser(0.55, 1.45)
    pad2.Update()
    line = ROOT.TLine(hr.GetBinCenter(1) - 0.5 * hr.GetBinWidth(1), 1, hr.GetBinCenter(hr.GetNbinsX()) + 0.5 * hr.GetBinWidth(hr.GetNbinsX()), 1.0)
    line.Draw()

    c.Print(os.path.join(options.input, "comparison_plots", f"{comp_name}_comparison.pdf"))
    c.Print(os.path.join(options.input, "comparison_plots", f"{comp_name}_comparison.png"))

    f.Close()


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-i', '--input',
                      action="store", dest="input",
                      default="", help="Input histogram folder")
    parser.add_option('-1', '--histogram1',
                      action="store", dest="h1",
                      default="", help="Sample histogram 1 to run over")
    parser.add_option('-2', '--histogram2',
                      action="store", dest="h2",
                      default="", help="Sample histogram 2 to run over")

    # parse input arguments
    options, args = parser.parse_args()

    if not os.path.isdir(os.path.join(options.input, "comparison_plots")):
        os.makedirs(os.path.join(options.input, "comparison_plots"))

    # run
    main(options, args)
