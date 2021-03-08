#!/usr/bin/env python
import os
import ROOT

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.SetAtlasStyle()

# Include custom PDFs
ROOT.gROOT.LoadMacro(os.path.join(os.path.dirname(__file__), "RooTwoSidedCBShape.cc"))

# Get the histogram
f = ROOT.TFile("wplusd_spg_comparison_norw/histograms_OS_0tag_Dplus.root", "READ")

# Make output folder
if not os.path.isdir("fits"):
    os.makedirs("fits")

for name in ['Wjets_emu_Matched_OS-SS_OS_0tag_Dplus_Dmeson_m', 'ForcedDecay_DPlus_Matched_OS_0tag_Dplus_Dmeson_m']:

    # print
    print(f"\n\nXXX {name} XXX\n\n")

    # Get the histogram
    h = f.Get(name)

    # Observable, i.e. the invariant mass
    x = ROOT.RooRealVar("x", "m(D^{+})", 1.7, 2.15, "GeV")
    x_arg_list = ROOT.RooArgList(x)
    x_arg_set = ROOT.RooArgSet(x)

    # Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
    dh = ROOT.RooDataHist("dh", "dh", x_arg_list, ROOT.RooFit.Import(h))

    # Mean, Sigma parameters
    mean = ROOT.RooRealVar("mean", "mean", 1.869, 1.869, 1.869)
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.02, 0.0, 1.0)

    # Signal
    nLo = ROOT.RooRealVar("nLo", "nLo", 1, 0, 50)
    nHi = ROOT.RooRealVar("nHi", "nHi", 1, 0, 50)
    aLo = ROOT.RooRealVar("aLo", "aLo", 1, 0, 50)
    aHi = ROOT.RooRealVar("aHi", "aHi", 1, 0, 50)
    sig = ROOT.RooTwoSidedCBShape("signal", "signal", x, mean, sigma, aLo, nLo, aHi, nHi)

    # Do the fit
    result = sig.fitTo(dh, ROOT.RooFit.Minimizer("Minuit2", "Migrad"),
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
    c = ROOT.TCanvas("c", "c", 800, 600)
    frame.Draw()
    c.Print(f"fits/{name}.pdf")
    frame.SetMaximum(1)
    frame.SetMinimum(1e-5)
    c.SetLogy()
    c.Print(f"fits/{name}_LOG.pdf")
