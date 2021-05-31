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

# Get the histogram
f = ROOT.TFile("/global/cscratch1/sd/mmuskinj/charmpp/v8/differential_fit_4_spg/wplusd_spg_comparison/histograms.root", "READ")

# Make output folder
if not os.path.isdir("fits"):
    os.makedirs("fits")

# histogram names
names = [
    "SPG_Matched_truth_pt_bin1_OS_Dplus_Matched_truth_pt_bin1_pt_bin1_Dmeson_m",
    "SPG_Matched_truth_pt_bin2_OS_Dplus_Matched_truth_pt_bin2_pt_bin2_Dmeson_m",
    "SPG_Matched_truth_pt_bin3_OS_Dplus_Matched_truth_pt_bin3_pt_bin3_Dmeson_m",
    "SPG_Matched_truth_pt_bin4_OS_Dplus_Matched_truth_pt_bin4_pt_bin4_Dmeson_m",
    "SPG_Matched_truth_pt_bin5_OS_Dplus_Matched_truth_pt_bin5_pt_bin5_Dmeson_m",
    "SPG_Matched_truth_pt_bin1_SS_Dplus_Matched_truth_pt_bin1_pt_bin1_Dmeson_m",
    "SPG_Matched_truth_pt_bin2_SS_Dplus_Matched_truth_pt_bin2_pt_bin2_Dmeson_m",
    "SPG_Matched_truth_pt_bin3_SS_Dplus_Matched_truth_pt_bin3_pt_bin3_Dmeson_m",
    "SPG_Matched_truth_pt_bin4_SS_Dplus_Matched_truth_pt_bin4_pt_bin4_Dmeson_m",
    "SPG_Matched_truth_pt_bin5_SS_Dplus_Matched_truth_pt_bin5_pt_bin5_Dmeson_m",
    "Wjets_emu_Matched_truth_pt_bin1_OS_Dplus_Matched_truth_pt_bin1_pt_bin1_Dmeson_m",
    "Wjets_emu_Matched_truth_pt_bin2_OS_Dplus_Matched_truth_pt_bin2_pt_bin2_Dmeson_m",
    "Wjets_emu_Matched_truth_pt_bin3_OS_Dplus_Matched_truth_pt_bin3_pt_bin3_Dmeson_m",
    "Wjets_emu_Matched_truth_pt_bin4_OS_Dplus_Matched_truth_pt_bin4_pt_bin4_Dmeson_m",
    "Wjets_emu_Matched_truth_pt_bin5_OS_Dplus_Matched_truth_pt_bin5_pt_bin5_Dmeson_m",
    "Wjets_emu_Matched_truth_pt_bin1_SS_Dplus_Matched_truth_pt_bin1_pt_bin1_Dmeson_m",
    "Wjets_emu_Matched_truth_pt_bin2_SS_Dplus_Matched_truth_pt_bin2_pt_bin2_Dmeson_m",
    "Wjets_emu_Matched_truth_pt_bin3_SS_Dplus_Matched_truth_pt_bin3_pt_bin3_Dmeson_m",
    "Wjets_emu_Matched_truth_pt_bin4_SS_Dplus_Matched_truth_pt_bin4_pt_bin4_Dmeson_m",
    "Wjets_emu_Matched_truth_pt_bin5_SS_Dplus_Matched_truth_pt_bin5_pt_bin5_Dmeson_m",
]


def main(options, args):

    # input file
    f = ROOT.TFile(options.input, "READ")

    # output file
    out = ROOT.TFile("fits/Dplus_mass_fit.root", "RECREATE")

    # loop
    for name in names:

        # print
        print(f"\n\nXXX {name} XXX\n\n")

        # Get the histogram
        f.cd()
        h_tmp = f.Get(name)
        h = h_tmp.Clone(f"{h_tmp.GetName()}_clone")
        h.Scale(100. / h.GetSumOfWeights())

        # Observable, i.e. the invariant mass
        x = ROOT.RooRealVar("x", "m(D^{+})", 1.7, 2.15, "GeV")
        x_arg_list = ROOT.RooArgList(x)

        # Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
        dh = ROOT.RooDataHist("dh", "dh", x_arg_list, ROOT.RooFit.Import(h))

        # Signal
        if "pt_bin5" not in name:
            mean = ROOT.RooRealVar("mean", "mean", 1.870, 1.870, 1.870)
            aLo = ROOT.RooRealVar("aLo", "aLo", 1.4, 1.4, 1.4)
            aHi = ROOT.RooRealVar("aHi", "aHi", 1.4, 1.4, 1.4)
        else:
            mean = ROOT.RooRealVar("mean", "mean", 1.873, 1.873, 1.873)
            aLo = ROOT.RooRealVar("aLo", "aLo", 1.35, 1.3, 1.3)
            aHi = ROOT.RooRealVar("aHi", "aHi", 1.35, 1.3, 1.3)
        sigma = ROOT.RooRealVar("sigma", "sigma", 0.02, 0.0, 0.05)
        nLo = ROOT.RooRealVar("nLo", "nLo", 1, 0, 100)
        nHi = ROOT.RooRealVar("nHi", "nHi", 1, 0, 100)
        sig = ROOT.RooTwoSidedCBShape("signal", "signal", x, mean, sigma, aLo, nLo, aHi, nHi)

        # Do the fit
        _ = sig.fitTo(dh, ROOT.RooFit.Minimizer("Minuit2", "Migrad"),
                      ROOT.RooFit.SumW2Error(ROOT.kTRUE),
                      ROOT.RooFit.Save(),
                      ROOT.RooFit.PrintLevel(1))

        s_component = ROOT.RooArgSet(sig)
        frame = x.frame(ROOT.RooFit.Title("Fit D+ mass"))

        # data
        dh.plotOn(frame, ROOT.RooFit.MarkerColor(ROOT.kBlack), ROOT.RooFit.Name("data"), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))

        # signal
        sig.plotOn(frame, ROOT.RooFit.Components(s_component), ROOT.RooFit.LineColor(
            ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name("s"))

        # canvas
        c = ROOT.TCanvas(name, name, 1000, 1000)
        frame.Draw()
        frame.SetMaximum(1.5 * h.GetMaximum())

        # pint text
        ROOT.ATLASLabel(0.17, 0.90, "Internal", 1)
        ROOT.myText(0.17, 0.9 - 1 * 0.06, 1, "#sqrt{s} = 13 TeV")
        ROOT.myText(0.17, 0.9 - 2 * 0.06, 1, f"D#rightarrowK#pi#pi {name.split('Matched_')[0]}")
        ROOT.myText(0.17, 0.9 - 3 * 0.06, 1, (name.split("_Dmeson_m")[0]).split("Matched_")[-1])
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
        frame.SetMaximum(1000 * h.GetMaximum())
        frame.SetMinimum(1e-3)
        c.SetLogy()
        c.Print(f"fits/{name}_LOG.pdf")

    f.Close()
    out.Close()


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
