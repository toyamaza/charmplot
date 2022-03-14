#!/usr/bin/env python
import os
import ROOT

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasLabels.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasUtils.C"))
ROOT.SetAtlasStyle()

# Include custom PDFs
ROOT.gROOT.LoadMacro(os.path.join(os.path.dirname(__file__), "RooTwoSidedCBShape.cc"))

# Make output folder
if not os.path.isdir("fits"):
    os.makedirs("fits")

# histogram names
names = [
    "MGPy8EG_NLO_WplusD_Matched_OS_0tag_Dstar_Dmeson_mdiff_norebin",
]


def main(options, args):

    # input file
    f = ROOT.TFile(options.input, "READ")

    # output file
    out = ROOT.TFile("fits/Dstar_mass_fit.root", "RECREATE")

    # output file for histograms
    out_histo = ROOT.TFile("DstarMassFit.root", "RECREATE")

    # loop
    for name in names:

        # print
        print(f"\n\nXXX {name} XXX\n\n")

        # Get the histogram
        f.cd()
        h_tmp = f.Get(name)
        h = h_tmp.Clone(f"{h_tmp.GetName()}_clone")
        h.Scale(100. / h.GetSumOfWeights())
        # h.Rebin(4)

        # Observable, i.e. the invariant mass
        x = ROOT.RooRealVar("x", "m(D*^{+}-D^{0})", 140, 180, "MeV")
        x_arg_list = ROOT.RooArgList(x)

        # Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
        dh = ROOT.RooDataHist("dh", "dh", x_arg_list, ROOT.RooFit.Import(h))

        # Signal
        mean = ROOT.RooRealVar("mean", "mean", 145.5, 144.0, 147.0)
        aLo = ROOT.RooRealVar("aLo", "aLo", 1, 0, 10)
        aHi = ROOT.RooRealVar("aHi", "aHi", 1, 0, 10)
        sigma = ROOT.RooRealVar("sigma", "sigma", 0.65, 0., 2.0)
        nLo = ROOT.RooRealVar("nLo", "nLo", 99, 0, 100)
        nHi = ROOT.RooRealVar("nHi", "nHi", 99, 0, 100)
        sig = ROOT.RooTwoSidedCBShape("signal", "signal", x, mean, sigma, aLo, nLo, aHi, nHi)

        # Do the fit
        _ = sig.fitTo(dh, ROOT.RooFit.Minimizer("Minuit2", "Migrad"),
                      ROOT.RooFit.SumW2Error(ROOT.kTRUE),
                      ROOT.RooFit.Save(),
                      ROOT.RooFit.PrintLevel(1))

        # up / down resolution
        sigma_up_val = sigma.getVal() * 1.05
        sigma_dn_val = sigma.getVal() * 0.95
        sigma_up = ROOT.RooRealVar("sigma_up", "sigma_up", sigma_up_val, sigma_up_val, sigma_up_val)
        sigma_dn = ROOT.RooRealVar("sigma_dn", "sigma_dn", sigma_dn_val, sigma_dn_val, sigma_dn_val)
        sig_up = ROOT.RooTwoSidedCBShape("signal_up", "signal_up", x, mean, sigma_up, aLo, nLo, aHi, nHi)
        sig_dn = ROOT.RooTwoSidedCBShape("signal_dn", "signal_dn", x, mean, sigma_dn, aLo, nLo, aHi, nHi)
        s_component = ROOT.RooArgSet(sig)
        frame = x.frame(ROOT.RooFit.Title("Fit D* mass"))

        # data
        dh.plotOn(frame, ROOT.RooFit.MarkerColor(ROOT.kBlack), ROOT.RooFit.Name("data"), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))

        # signal
        sig.plotOn(frame, ROOT.RooFit.Components(s_component), ROOT.RooFit.LineColor(
            ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name("s"))

        # signal up/dn
        sig_up.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name("s_up"))
        sig_dn.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name("s_dn"))

        # canvas
        c = ROOT.TCanvas(name, name, 1000, 1000)
        frame.Draw()
        frame.SetMaximum(1.5 * h.GetMaximum())

        # print text
        ROOT.ATLASLabel(0.17, 0.90, "Internal", 1)
        ROOT.myText(0.17, 0.9 - 1 * 0.06, 1, "#sqrt{s} = 13 TeV")
        ROOT.myText(0.17, 0.9 - 2 * 0.06, 1, f"D*#rightarrowK#pi#pi {name.split('_WplusD_Matched_')[0]}")
        ROOT.myText(0.17, 0.9 - 3 * 0.06, 1, (name.split("_Dmeson_mdiff")[0]).split("_WplusD_Matched_")[-1])
        ROOT.myText(0.68, 0.9 - 0 * 0.06, 1, f"aHi: {aHi.getVal():.4f}")
        ROOT.myText(0.68, 0.9 - 1 * 0.06, 1, f"aLo: {aLo.getVal():.4f}")
        ROOT.myText(0.68, 0.9 - 2 * 0.06, 1, f"mean: {mean.getVal():.4f}")
        ROOT.myText(0.68, 0.9 - 3 * 0.06, 1, f"nHi: {nHi.getVal():.4f}")
        ROOT.myText(0.68, 0.9 - 4 * 0.06, 1, f"nLo: {nLo.getVal():.4f}")
        ROOT.myText(0.68, 0.9 - 5 * 0.06, 1, f"sigma: {sigma.getVal():.4f}")

        # save to output
        h_out = ROOT.TH1D(f"pars_{name}", f"pars_{name}", 6, -0.5, 5.5)
        h_out.SetBinContent(1, aHi.getVal())
        h_out.SetBinContent(2, aLo.getVal())
        h_out.SetBinContent(3, mean.getVal())
        h_out.SetBinContent(4, nHi.getVal())
        h_out.SetBinContent(5, nLo.getVal())
        h_out.SetBinContent(6, sigma.getVal())
        h_out.SetBinError(1, aHi.getError())
        h_out.SetBinError(2, aLo.getError())
        h_out.SetBinError(3, mean.getError())
        h_out.SetBinError(4, nHi.getError())
        h_out.SetBinError(5, nLo.getError())
        h_out.SetBinError(6, sigma.getError())
        h_out.GetXaxis().SetBinLabel(1, aHi.GetName())
        h_out.GetXaxis().SetBinLabel(2, aLo.GetName())
        h_out.GetXaxis().SetBinLabel(3, mean.GetName())
        h_out.GetXaxis().SetBinLabel(4, nHi.GetName())
        h_out.GetXaxis().SetBinLabel(5, nLo.GetName())
        h_out.GetXaxis().SetBinLabel(6, sigma.GetName())
        out.cd()
        h_out.Write()

        c.Print(f"fits/{name}.pdf")
        c.Print(f"fits/{name}.png")
        # frame.SetMaximum(1000 * h.GetMaximum())
        # frame.SetMinimum(1e-3)
        # c.SetLogy()
        # c.Print(f"fits/{name}_LOG.pdf")
        # c.Print(f"fits/{name}_LOG.png")

        # save mass histograms
        out_histo.cd()
        f_sig = sig.asTF(x)
        f_sig_up = sig_up.asTF(x)
        f_sig_dn = sig_dn.asTF(x)
        h_sig = h.Clone(name.replace("_Dmeson_mdiff_norebin", "__Dmeson_mdiff"))
        h_sig_up = h.Clone("sigma_up_" + name.replace("_Dmeson_mdiff_norebin", "__Dmeson_mdiff"))
        h_sig_dn = h.Clone("sigma_dn_" + name.replace("_Dmeson_mdiff_norebin", "__Dmeson_mdiff"))
        for i in range(1, h_sig.GetNbinsX() + 1):
            h_sig.SetBinContent(i, f_sig.Integral(h.GetBinCenter(i) - 0.5 * h.GetBinWidth(i), h.GetBinCenter(i) + 0.5 * h.GetBinWidth(i)) / h.GetBinWidth(i))
            h_sig_up.SetBinContent(i, f_sig_up.Integral(h.GetBinCenter(i) - 0.5 * h.GetBinWidth(i), h.GetBinCenter(i) + 0.5 * h.GetBinWidth(i)) / h.GetBinWidth(i))
            h_sig_dn.SetBinContent(i, f_sig_dn.Integral(h.GetBinCenter(i) - 0.5 * h.GetBinWidth(i), h.GetBinCenter(i) + 0.5 * h.GetBinWidth(i)) / h.GetBinWidth(i))
        h_sig.Write()
        h_sig_up.Write()
        h_sig_dn.Write()

    f.Close()
    out.Close()
    out_histo.Close()


if __name__ == "__main__":


    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-i', '--input',
                      action="store", dest="input",
                      help="Path to the histograms.root file from the initial spg comparison.")

    # parse input arguments
    options, args = parser.parse_args()

    # run
    main(options, args)
