#!/usr/bin/env python
import os
import ROOT
import re

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasLabels.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasUtils.C"))
ROOT.SetAtlasStyle()

# Include custom PDFs
# ROOT.gROOT.LoadMacro(os.path.join(os.path.dirname(__file__), "RooTwoSidedCBShape.cc"))

# Make output folder
if not os.path.isdir("fits_bkg"):
    os.makedirs("fits_bkg")

# histogram names
names = [
    "Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_Dmeson_mdiff",
    "Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin1_Dmeson_mdiff",
    "Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin2_Dmeson_mdiff",
    "Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin3_Dmeson_mdiff",
    "Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin4_Dmeson_mdiff",
    "Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin5_Dmeson_mdiff",
    "Wjets_emu_Rest_PostProc_Loose_OS_0tag_Dstar_Wjets_emu_Rest_Dmeson_mdiff",
    "Wjets_emu_Rest_PostProc_Loose_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin1_Dmeson_mdiff",
    "Wjets_emu_Rest_PostProc_Loose_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin2_Dmeson_mdiff",
    "Wjets_emu_Rest_PostProc_Loose_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin3_Dmeson_mdiff",
    "Wjets_emu_Rest_PostProc_Loose_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin4_Dmeson_mdiff",
    "Wjets_emu_Rest_PostProc_Loose_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin5_Dmeson_mdiff",
    # "Sherpa_Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_Dmeson_mdiff",
    # "Sherpa_Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin1_Dmeson_mdiff",
    # "Sherpa_Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin2_Dmeson_mdiff",
    # "Sherpa_Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin3_Dmeson_mdiff",
    # "Sherpa_Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin4_Dmeson_mdiff",
    # "Sherpa_Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin5_Dmeson_mdiff",
    # "Powheg_Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_Dmeson_mdiff",
    # "Powheg_Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin1_Dmeson_mdiff",
    # "Powheg_Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin2_Dmeson_mdiff",
    # "Powheg_Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin3_Dmeson_mdiff",
    # "Powheg_Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin4_Dmeson_mdiff",
    # "Powheg_Wjets_emu_Rest_PostProc_OS_0tag_Dstar_Wjets_emu_Rest_pt_bin5_Dmeson_mdiff",
]


def main(options, args):

    # input file
    f = ROOT.TFile(options.input, "READ")

    # output file
    out = ROOT.TFile("fits_bkg/Dstar_mass_fit_bkg.root", "RECREATE")
    out_gen = ROOT.TFile("wjets_bkg_fit.root", "RECREATE")
    out_var = ROOT.TFile("wjets_bkg_fit_var.root", "RECREATE")
    out_sherpa = ROOT.TFile("sherpa_wjets_bkg_fit.root", "RECREATE")
    out_powheg = ROOT.TFile("powheg_wjets_bkg_fit.root", "RECREATE")

    # Create pdf
    ROOT.RooClassFactory.makePdf("bkg_pdf", "x,a0,a1,a2", "", "a0+a1*log(x-a2)")

    # Compile pdf class
    ROOT.gROOT.ProcessLineSync(".x bkg_pdf.cxx+")

    # loop
    for name in names:

        # print
        print(f"\n\nXXX {name} XXX\n\n")

        # Get pt bin
        if re.search("pt_bin", name):
            pt_bin = (name.split("_Dmeson")[0]).split("Rest_")[-1]
        else:
            pt_bin = ''

        # Get pt bin
        if re.search("pt_bin", name):
            pt_bin = (name.split("_Dmeson")[0]).split("Rest_")[-1]
        else:
            pt_bin = ''

        # Get Sherpa or Madgraph
        if re.search("Sherpa", name):
            sherpa = True
            powheg = False
            mc_string = "Sherpa_"
        elif re.search("Powheg", name):
            powheg = True
            sherpa = False
            mc_string = "Powheg_"
        else:
            sherpa = False
            powheg = False

        if re.search("Loose", name):
            loose = True
        else:
            loose = False

        # Get the histogram
        f.cd()
        h_tmp = f.Get(name)
        h = h_tmp.Clone(f"{h_tmp.GetName()}_clone")
        # h.Rebin(10)
        # h.Scale(100. / h.GetSumOfWeights())

        # Observable, i.e. the invariant mass
        x = ROOT.RooRealVar("x", "m(D*^{+}-D^{0})", 140, 180, "MeV")
        x_arg_list = ROOT.RooArgList(x)

        # Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
        dh = ROOT.RooDataHist("dh", "dh", x_arg_list, ROOT.RooFit.Import(h))

        # Signal
        a0 = ROOT.RooRealVar("a0", "", 0., -100., 100.)
        a1 = ROOT.RooRealVar("a1", "", 1., -100., 100.)
        a2 = ROOT.RooRealVar("a2", "", 135., 130., 140.)
        bkg_fn = ROOT.bkg_pdf("bkg_fn", "bkg_fn", x, a0, a1, a2)

        # Do the fit
        fit_result = bkg_fn.fitTo(dh, ROOT.RooFit.Minimizer("Minuit2", "Migrad"),
                                  ROOT.RooFit.SumW2Error(ROOT.kTRUE),
                                  ROOT.RooFit.Save(),
                                  ROOT.RooFit.PrintLevel(1))

        s_component = ROOT.RooArgSet(bkg_fn)
        frame = x.frame(ROOT.RooFit.Title("Fit D* mass"))

        # data
        dh.plotOn(frame, ROOT.RooFit.MarkerColor(ROOT.kBlack), ROOT.RooFit.Name("data"), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))

        # signal
        bkg_fn.plotOn(frame, ROOT.RooFit.Components(s_component), ROOT.RooFit.LineColor(
            ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name("s"))

        # Extract covariance and correlation matrix as ROOT.TMatrixDSym
        cor = fit_result.correlationMatrix()
        cov = fit_result.covarianceMatrix()
        # Print correlation, matrix
        print("correlation matrix")
        cor.Print()
        print("covariance matrix")
        cov.Print()

        # canvas
        c = ROOT.TCanvas(name, name, 1000, 1000)
        frame.Draw()
        frame.SetMaximum(1.5 * h.GetMaximum())

        # pint text
        ROOT.ATLASLabel(0.17, 0.90, "Internal", 1)
        ROOT.myText(0.17, 0.9 - 1 * 0.06, 1, "#sqrt{s} = 13 TeV")
        ROOT.myText(0.17, 0.9 - 2 * 0.06, 1, "D*#rightarrowK#pi#pi")
        ROOT.myText(0.17, 0.9 - 3 * 0.06, 1, (name.split("_Wjets")[0]).split("PostProc_")[-1])
        ROOT.myText(0.17, 0.9 - 4 * 0.06, 1, (name.split("_Dmeson_mdiff")[0]).split("Dstar_")[-1])
        ROOT.myText(0.68, 0.9 - 0 * 0.06, 1, f"a0: {a0.getVal():.4f}")
        ROOT.myText(0.68, 0.9 - 1 * 0.06, 1, f"a1: {a1.getVal():.4f}")
        ROOT.myText(0.68, 0.9 - 2 * 0.06, 1, f"a2: {a2.getVal():.4f}")

        # save to output
        h_out = ROOT.TH1D(f"pars_{name}", f"pars_{name}", 3, -0.5, 5.5)
        h_out.SetBinContent(1, a0.getVal())
        h_out.SetBinContent(2, a1.getVal())
        h_out.SetBinContent(3, a2.getVal())
        h_out.SetBinError(1, a0.getError())
        h_out.SetBinError(2, a1.getError())
        h_out.SetBinError(3, a2.getError())
        h_out.GetXaxis().SetBinLabel(1, a0.GetName())
        h_out.GetXaxis().SetBinLabel(2, a1.GetName())
        h_out.GetXaxis().SetBinLabel(3, a2.GetName())
        out.cd()
        h_out.Write()

        if pt_bin:
            c.Print(f"fits_bkg/{name}_{pt_bin}.pdf")
            c.Print(f"fits_bkg/{name}_{pt_bin}.png")
        else:
            c.Print(f"fits_bkg/{name}.pdf")
            c.Print(f"fits_bkg/{name}.png")

        # frame.SetMaximum(1000 * h.GetMaximum())
        # frame.SetMinimum(1e-3)
        # c.SetLogy()
        # c.Print(f"fits/{name}_LOG.pdf")
        # c.Print(f"fits/{name}_LOG.png")

        # Convert fit to TF1
        wjets_fn = ROOT.TF1("wjets_fn", "[0] + [1]*log(x-[2])", 140, 180)
        wjets_fn.SetParameters(a0.getVal(), a1.getVal(), a2.getVal())
        wjets_fn.SetParNames("a0", "a1", "a2")
        wjets_fn.SetParError(0, a0.getError())
        wjets_fn.SetParError(1, a1.getError())
        wjets_fn.SetParError(2, a2.getError())

        # Fit to input histogram
        fit_result_wjets = h.Fit("wjets_fn", "S")

        # Create histogram generated from Fit
        # h_gen = ROOT.TH1F("h_gen", "h_gen", h.GetNbinsX(),140,180)
        h_gen = h.Clone("h_gen")

        # Fill histogram
        for i in range(1, h_gen.GetNbinsX() + 1):
            # Do not remove this print statement (for whatever reason it makes everything work, its...magic)
            print(f'BinIntegralError = {wjets_fn.IntegralError(h.GetBinCenter(i) - 0.5 * h.GetBinWidth(i), h.GetBinCenter(i) + 0.5 * h.GetBinWidth(i), wjets_fn.GetParameters(), fit_result_wjets.GetCovarianceMatrix().GetMatrixArray())}')
            h_gen.SetBinContent(i, wjets_fn.Integral(h.GetBinCenter(i) - 0.5 * h.GetBinWidth(i), h.GetBinCenter(i) + 0.5 * h.GetBinWidth(i)) / h.GetBinWidth(i))
            h_gen.SetBinError(i, wjets_fn.IntegralError(h.GetBinCenter(i) - 0.5 * h.GetBinWidth(i), h.GetBinCenter(i) + 0.5 * h.GetBinWidth(i), wjets_fn.GetParameters(), fit_result_wjets.GetCovarianceMatrix().GetMatrixArray()) / h.GetBinWidth(i))

        # Create the variational histograms
        # h_var_up = ROOT.TH1F("h_var_up", "h_var_up", h_gen.GetNbinsX(),140,180)
        h_var_up = h_gen.Clone("h_var_up")

        # Fill variation histogram
        for i in range(1, h_gen.GetNbinsX() + 1):
            h_var_up.SetBinContent(i, h_gen.GetBinContent(i) + h_gen.GetBinError(i))

        # Get Confidence Intervals
        print(fit_result_wjets.GetCovarianceMatrix().Print())
        hConf1s = h_gen.Clone("hConf1s")
        hConf2s = h_gen.Clone("hConf2s")
        hConf1s.Reset()
        hConf2s.Reset()
        ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(hConf1s, 0.68)
        ROOT.TVirtualFitter.GetFitter().GetConfidenceIntervals(hConf2s, 0.95)

        # Configure Confidence Interval Colors
        hConf1s.SetStats(0)
        hConf2s.SetStats(0)
        hConf1s.SetMarkerSize(0)
        hConf2s.SetMarkerSize(0)
        hConf1s.SetFillColor(ROOT.kGreen + 1)
        hConf2s.SetFillColor(ROOT.kOrange)
        hConf1s.SetLineColor(ROOT.kGreen + 1)
        hConf2s.SetLineColor(ROOT.kOrange)

        # Plot generated histogram from fit
        c_gen = ROOT.TCanvas(name + '_gen', name + '_gen', 1000, 1000)
        h_gen.SetMaximum(1.5 * h_gen.GetMaximum())
        h_gen.Draw("E0")
        hConf2s.Draw("e3 same")
        hConf1s.Draw("e3 same")
        h.Draw("P* same")
        ROOT.ATLASLabel(0.17, 0.90, "Internal", 1)
        ROOT.myText(0.17, 0.9 - 1 * 0.06, 1, "#sqrt{s} = 13 TeV")
        ROOT.myText(0.17, 0.9 - 2 * 0.06, 1, "D*#rightarrowK#pi#pi Fit")
        ROOT.myText(0.17, 0.9 - 3 * 0.06, 1, (name.split("_Wjets")[0]).split("PostProc_")[-1])
        ROOT.myText(0.17, 0.9 - 4 * 0.06, 1, (name.split("_Dmeson_mdiff")[0]).split("Dstar_")[-1])
        if sherpa:
            out_sherpa.cd()
        elif powheg:
            out_powheg.cd()
        else:
            out_gen.cd()

        if sherpa or powheg:
            if pt_bin:
                h_gen.Write(f"{mc_string}Wjets_emu_Rest_Fit_{pt_bin}__Dmeson_mdiff")
            else:
                h_gen.Write(f"{mc_string}Wjets_emu_Rest_Fit__Dmeson_mdiff")
        elif loose:
            if pt_bin:
                h_gen.Write(f"Wjets_emu_Rest_Fit_Loose_{pt_bin}__Dmeson_mdiff")
            else:
                h_gen.Write("Wjets_emu_Rest_Fit_Loose__Dmeson_mdiff")
        else:
            if pt_bin:
                h_gen.Write(f"Wjets_emu_Rest_Fit_{pt_bin}__Dmeson_mdiff")
            else:
                h_gen.Write("Wjets_emu_Rest_Fit__Dmeson_mdiff")

        if pt_bin:
            c_gen.Print(f"fits_bkg/{name}_{pt_bin}_gen.pdf")
            c_gen.Print(f"fits_bkg/{name}_{pt_bin}_gen.png")
        else:
            c_gen.Print(f"fits_bkg/{name}_{pt_bin}_gen.pdf")
            c_gen.Print(f"fits_bkg/{name}_{pt_bin}_gen.png")

        # Plot upward variation
        c_var_up = ROOT.TCanvas(name + '_1up', name + '_1up', 1000, 1000)
        h_var_up.SetMaximum(1.5 * h_var_up.GetMaximum())
        h_var_up.Draw("hist")
        ROOT.ATLASLabel(0.17, 0.90, "Internal", 1)
        ROOT.myText(0.17, 0.9 - 1 * 0.06, 1, "#sqrt{s} = 13 TeV")
        ROOT.myText(0.17, 0.9 - 2 * 0.06, 1, "D*#rightarrowK#pi#pi 1up")
        ROOT.myText(0.17, 0.9 - 3 * 0.06, 1, (name.split("_Wjets")[0]).split("PostProc_")[-1])
        ROOT.myText(0.17, 0.9 - 4 * 0.06, 1, (name.split("_Dmeson_mdiff")[0]).split("Dstar_")[-1])
        out_var.cd()
        if sherpa or powheg:
            if pt_bin:
                h_var_up.Write(f"{mc_string}Wjets_emu_Rest_Fit_{pt_bin}__Dmeson_mdiff")
            else:
                h_var_up.Write(f"{mc_string}Wjets_emu_Rest_Fit__Dmeson_mdiff")
        elif not loose:
            if pt_bin:
                h_var_up.Write(f"Wjets_emu_Rest_Fit_{pt_bin}__Dmeson_mdiff")
            else:
                h_var_up.Write("Wjets_emu_Rest_Fit__Dmeson_mdiff")

        if pt_bin:
            c_var_up.Print(f"fits_bkg/{name}_{pt_bin}_1up.pdf")
            c_var_up.Print(f"fits_bkg/{name}_{pt_bin}_1up.png")
        else:
            c_var_up.Print(f"fits_bkg/{name}_{pt_bin}_1up.pdf")
            c_var_up.Print(f"fits_bkg/{name}_{pt_bin}_1up.png")

    out_var.Close()
    out_gen.Close()
    out_sherpa.Close()
    out_powheg.Close()
    out.Close()
    f.Close()


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
