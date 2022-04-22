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


def createCanvasPads(name):
    c = ROOT.TCanvas(name, name, 1200, 1200)
    # Upper histogram plot is pad1
    pad1 = ROOT.TPad(f"pad1_{name}", f"pad1_{name}", 0, 0.30, 1, 1.0)
    pad1.SetTopMargin(0.1)
    pad1.SetBottomMargin(0.05 / 0.70)
    pad1.SetFillStyle(4000)
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()
    pad2 = ROOT.TPad(f"pad2_{name}", f"pad2_{name}", 0, 0.00, 1, 0.40)
    pad2.SetTopMargin(0.05 / 0.40)
    pad2.SetBottomMargin(0.10 / 0.40)
    pad2.SetFillStyle(4000)
    pad2.Draw()

    return c, pad1, pad2


def main(options):
    f = ROOT.TFile(f"{options.input}/outpoot.root", "READ")

    for pt_bin in range(1, 6):

        h_minus = f.Get(f"h_tot_postFit_OS_lep_minus_0tag_Dplus_pt_bin{pt_bin}_SS_lep_minus_0tag_Dplus_pt_bin{pt_bin}")
        h_plus = f.Get(f"h_tot_postFit_OS_lep_plus_0tag_Dplus_pt_bin{pt_bin}_SS_lep_plus_0tag_Dplus_pt_bin{pt_bin}")

        gr_minus = f.Get(f"h_tot_postFit_OS_lep_minus_0tag_Dplus_pt_bin{pt_bin}_SS_lep_minus_0tag_Dplus_pt_bin{pt_bin}_data")
        gr_plus = f.Get(f"h_tot_postFit_OS_lep_plus_0tag_Dplus_pt_bin{pt_bin}_SS_lep_plus_0tag_Dplus_pt_bin{pt_bin}_data")

        gr_minus.SetMarkerColor(ROOT.kRed)
        gr_minus.SetLineColor(ROOT.kRed)
        gr_plus.SetMarkerColor(ROOT.kBlue)
        gr_plus.SetLineColor(ROOT.kBlue)

        gr_ratio = gr_minus.Clone()
        for i in range(gr_ratio.GetN()):
            y_minus = gr_minus.GetY()[i] / h_minus.Integral()
            y_plus = gr_plus.GetY()[i] / h_plus.Integral()
            y_err_minus = gr_minus.GetEY()[i] / h_minus.Integral()
            y_err_plus = gr_plus.GetEY()[i] / h_plus.Integral()
            gr_minus.GetY()[i] = y_minus
            gr_plus.GetY()[i] = y_plus
            gr_minus.GetEY()[i] = y_err_minus
            gr_plus.GetEY()[i] = y_err_plus
            gr_minus.GetEX()[i] = h_minus.GetBinWidth(i + 1) / 2.
            gr_plus.GetEX()[i] = h_minus.GetBinWidth(i + 1) / 2.

            y = y_plus / y_minus
            y_err = y * ((y_err_minus / y_minus)**2 + (y_err_plus / y_plus)**2)**(0.5)

            gr_ratio.SetPoint(i, gr_minus.GetX()[i], y)
            gr_ratio.SetPointError(i, h_minus.GetBinWidth(i + 1) / 2., y_err)

        # proxy histogram
        h_proxy_up = h_minus.Clone(f"{h_minus.GetName()}_proxy_up")
        h_proxy_dn = h_minus.Clone(f"{h_minus.GetName()}_proxy_dn")
        for i in range(1, h_proxy_up.GetNbinsX() + 1):
            h_proxy_up.SetBinContent(i, -1e5)
            h_proxy_up.SetBinError(i, 0)
            h_proxy_dn.SetBinContent(i, -1e5)
            h_proxy_dn.SetBinError(i, 0)

        GLOBAL_SF = 1.0
        h_proxy_up.GetYaxis().SetTitleSize(h_proxy_up.GetYaxis().GetTitleSize() * GLOBAL_SF)
        h_proxy_up.GetYaxis().SetTitleOffset(h_proxy_up.GetYaxis().GetTitleOffset() * (1 / (GLOBAL_SF * 1.1)))
        h_proxy_up.GetYaxis().SetLabelSize(h_proxy_up.GetYaxis().GetLabelSize() * GLOBAL_SF)
        h_proxy_up.GetXaxis().SetLabelSize(0)

        SF = 0.70 / 0.40
        h_proxy_dn.GetYaxis().SetTitleSize(h_proxy_dn.GetYaxis().GetTitleSize() * SF * GLOBAL_SF)
        h_proxy_dn.GetYaxis().SetTitleOffset(h_proxy_dn.GetYaxis().GetTitleOffset() * (1 / (GLOBAL_SF * 1.06 * SF)))
        h_proxy_dn.GetXaxis().SetTitleSize(h_proxy_dn.GetXaxis().GetTitleSize() * SF * GLOBAL_SF)
        h_proxy_dn.GetXaxis().SetTitleOffset(h_proxy_dn.GetXaxis().GetTitleOffset() * (1 / (GLOBAL_SF * 0.6 * SF)))
        h_proxy_dn.GetXaxis().SetLabelOffset(h_proxy_dn.GetXaxis().GetLabelOffset() * 4.0)
        h_proxy_dn.GetYaxis().SetLabelSize(h_proxy_dn.GetYaxis().GetLabelSize() * SF * GLOBAL_SF)
        h_proxy_dn.GetXaxis().SetLabelSize(h_proxy_dn.GetXaxis().GetLabelSize() * SF * GLOBAL_SF * 1.2)

        # tick marks
        h_proxy_up.GetYaxis().SetNdivisions(506)
        h_proxy_dn.GetYaxis().SetNdivisions(306)

        # title and range
        h_proxy_up.GetYaxis().SetRangeUser(1e-6, 0.5)
        h_proxy_up.GetYaxis().SetTitle("Normalized Entries / GeV")

        h_proxy_dn.GetYaxis().SetRangeUser(0.51, 1.49)
        h_proxy_dn.GetYaxis().SetTitle("Ratio W+ / W-")
        h_proxy_dn.GetXaxis().SetTitle("m(D#pm) [GeV]")

        # legend
        N = 2
        leg = ROOT.TLegend(0.58, 0.848 - N / 2. * (0.095), 0.92, 0.848)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(42)
        leg.SetTextFont(43)
        leg.AddEntry(gr_minus, "Data W-", "pe")
        leg.AddEntry(gr_plus, "Data W+", "pe")

        # tick marks
        gr_minus.GetYaxis().SetNdivisions(506)
        gr_ratio.GetYaxis().SetNdivisions(306)

        gr_ratio.GetYaxis().SetRangeUser(0.51, 1.49)
        gr_ratio.GetYaxis().SetTitle("Ratio W+ / W-")
        gr_ratio.GetXaxis().SetTitle("m(D) [GeV]")

        # multigraph
        mg = ROOT.TMultiGraph()
        mg.Add(gr_minus, "pe0")
        mg.Add(gr_plus, "pe0")

        # title
        mg.GetYaxis().SetTitle("Normalized Entries / GeV")

        c, pad1, pad2 = createCanvasPads(f"pt_bin{pt_bin}")
        pad1.cd()
        pad1.SetGridx()
        h_proxy_up.Draw()
        mg.Draw("pe")

        # ATLAS label
        l1 = ROOT.TLatex()
        l1.SetTextFont(73)
        l1.SetTextSize(42)
        l1.DrawLatexNDC(0.18, 0.82, "ATLAS")
        l2 = ROOT.TLatex()
        l2.SetTextFont(43)
        l2.SetTextSize(42)
        l2.DrawLatexNDC(0.33, 0.82, "Internal")
        l2.DrawLatexNDC(0.18, 0.82 - 1 * 0.05, "W#rightarrowl#nu+D, D#rightarrowK#pi#pi, OS-SS")
        l2.DrawLatexNDC(0.18, 0.82 - 2 * 0.05, f"pt_bin{pt_bin}")

        # legend
        leg.Draw()

        pad2.cd()
        pad2.SetGridy()
        pad2.SetGridx()
        h_proxy_dn.Draw()
        gr_ratio.Draw("pe")

        c.Print(f"{options.input}/pt_bin{pt_bin}_ratio.png")
        c.Print(f"{options.input}/pt_bin{pt_bin}_ratio.pdf")


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-i', '--input',
                      action="store", dest="input")

    # parse input arguments
    options, _ = parser.parse_args()

    # run
    main(options)
