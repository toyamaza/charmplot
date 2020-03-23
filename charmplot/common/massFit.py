from charmplot.common import canvas
from charmplot.common import utils
from charmplot.control import channel
from charmplot.control import variable
import logging
import os
import ROOT


# logging
logger = logging.getLogger(__name__)

# proxy histograms for legend
total_pdf_proxy = ROOT.TH1F("total_pdf_proxy", "", 1, 0, 1)
total_pdf_proxy.SetLineColor(ROOT.kBlue)
total_pdf_proxy.SetLineStyle(ROOT.kDashed)
background_proxy = ROOT.TH1F("background", "", 1, 0, 1)
background_proxy.SetLineColor(ROOT.kBlue)
signal_proxy = ROOT.TH1F("signal", "", 1, 0, 1)
signal_proxy.SetLineColor(ROOT.kRed)
mc_histograms = [
    total_pdf_proxy,
    background_proxy,
    signal_proxy,
]
mc_names = [
    "Total PDF",
    "Bkg. PDF",
    "Signal PDF",
]


class MassFit(object):

    def __init__(self, channel: channel.Channel, variable: variable.Variable, histogram: ROOT.TH1, label: str, output: str = ""):

        self.channel = channel
        self.name = channel.name
        self.fit_config = channel.mass_fit
        self.variable = variable
        self.h = histogram.Clone(f"{histogram.GetName()}_mass_fit")
        self.label = label
        self.output = output

        self.range = None
        self.bkg_function = None
        self.SumW2Error = None
        if "range" in self.fit_config.keys():
            self.range = self.fit_config["range"]
        if "bkg_function" in self.fit_config.keys():
            self.bkg_function = self.fit_config["bkg_function"]
        if "SumW2Error" in self.fit_config.keys():
            self.SumW2Error = self.fit_config["SumW2Error"]
        else:
            self.SumW2Error = True

        utils.set_to_positive(self.h)

    def fit(self):
        # Observable, i.e. the invariant mass
        x = ROOT.RooRealVar("x", "m(D^{+})", self.range[0], self.range[1], "GeV")
        x_arg_list = ROOT.RooArgList(x)

        # Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
        dh = ROOT.RooDataHist("dh", "dh", x_arg_list, ROOT.RooFit.Import(self.h))

        # Mean, Sigma parameters
        mean = ROOT.RooRealVar("mean", "mean", self.fit_config["mean"][0], self.fit_config["mean"][1], self.fit_config["mean"][2])
        sigma = ROOT.RooRealVar("sigma", "sigma", self.fit_config["sigma"][0], self.fit_config["sigma"][1], self.fit_config["sigma"][2])

        # Gaussian signal
        sig = ROOT.RooGaussian("signal", "signal", x, mean, sigma)

        if not self.bkg_function or self.bkg_function == "poly":
            # Polynomial for background
            a = []
            for i in range(self.fit_config["poly_order"]):
                ai = ROOT.RooRealVar(f"a{i}", f"a{i}", 1, -10, 10)
                a.append(ai)
            bkg = ROOT.RooPolynomial("bkg", "bkg", x, ROOT.RooArgList(*a), 0)
        elif self.bkg_function == "expo":
            # Exponential for background
            lambd = ROOT.RooRealVar("lambda", "slope", -0.1, -5., 0.)
            bkg = ROOT.RooExponential("bkg", "bkg", x, lambd)
        elif self.bkg_function == "epoly":
            args = [x]
            expr = []
            for i in range(self.fit_config["poly_order"]):
                ai = ROOT.RooRealVar(f"a{i+1}", f"a{i+1}", 1, -100, 100)
                args.append(ai)
                expr_i = f"a{i+1}"
                for j in range(i+1):
                    expr_i += " * x"
                expr += [expr_i]
            expr = " + ".join(expr)
            print (f"exp({expr})")
            bkg = ROOT.RooClassFactory.makePdfInstance(f"bkg_{self.channel}_{self.label}", f"exp({expr})", ROOT.RooArgList(*args))
        elif self.bkg_function == "power":
            args = [x]
            expr = []
            for i in range(self.fit_config["poly_order"]):
                ai = ROOT.RooRealVar(f"a{i}", f"a{i}", 1, -100, 100)
                args.append(ai)
                expr_i = f"a{i}"
                for j in range(i):
                    expr_i += " * x"
                expr += [expr_i]
            expr = " + ".join(expr)
            print (f"pow(x, {expr})")
            bkg = ROOT.RooClassFactory.makePdfInstance(f"bkg_{self.channel}_{self.label}", f"pow(x, {expr})", ROOT.RooArgList(*args))

        # Sum of two PDFs
        nsig = ROOT.RooRealVar("nsig", "nsig", 50000, 0, 1e6)
        nbkg = ROOT.RooRealVar("nbkg", "nbkg", 50000, 0, 1e6)
        total_pdf = ROOT.RooAddPdf("tot", "Total PDF", ROOT.RooArgList(sig, bkg), ROOT.RooArgList(nsig, nbkg))

        # Do the fit
        _ = total_pdf.fitTo(dh, ROOT.RooFit.Minimizer("Minuit2", "Migrad"),
                            ROOT.RooFit.SumW2Error(self.SumW2Error),
                            ROOT.RooFit.Save(),
                            ROOT.RooFit.PrintLevel(9))

        # Plot result
        s_component = ROOT.RooArgSet(sig)
        b_component = ROOT.RooArgSet(bkg)
        frame = x.frame()

        # data
        dh.plotOn(frame, ROOT.RooFit.MarkerColor(ROOT.kBlack), ROOT.RooFit.Name("data"), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))

        # background
        total_pdf.plotOn(frame, ROOT.RooFit.Components(b_component), ROOT.RooFit.LineColor(
            ROOT.kBlue), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name("b"))

        # total pdf
        total_pdf.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("total"))

        # signal
        total_pdf.plotOn(frame, ROOT.RooFit.Components(s_component), ROOT.RooFit.LineColor(
            ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name("s"))

        # Chi2
        floating = total_pdf.getParameters(dh).selectByAttrib("Constant", ROOT.kFALSE).getSize()
        ndf = dh.numEntries() - floating
        chi2 = frame.chiSquare("total", "data", floating)
        prob = ROOT.TMath.Prob(chi2 * ndf, ndf)

        # TCanvas object
        canv = canvas.CanvasMassFit(self.channel, self.variable, 800, 800, 0.0)
        canv.construct(self.h)

        # Draw frame
        canv.proxy_up.GetXaxis().SetRangeUser(self.range[0], self.range[1])
        canv.pad1.cd()
        frame.Draw("same")

        # Legend and additional text
        canv.text(f"{self.label} fit")
        canv.make_legend(self.h, mc_histograms=mc_histograms, mc_names=mc_names)
        canv.text_right("#chi^{2} = %.2f (%.1f%s)" % (chi2, prob * 100, "%"))
        canv.text_right("N_{S} = %.0f #pm %.0f" % (nsig.getVal(), nsig.getError()))
        canv.text_right("m = %.1f #pm %.2f" % (1000 * mean.getVal(), 1000 * mean.getError()))
        canv.text_right("#sigma = %.1f #pm %.2f" % (1000 * sigma.getVal(), 1000 * sigma.getError()))

        # set maximum after creating legend
        canv.set_maximum([self.h], self.variable)
        self.save_results(canv)

    def save_results(self, c):
        if not os.path.isdir(self.output):
            os.makedirs(self.output)
        c.canv.Print(os.path.join(self.output, f"{self.channel}_{self.variable.name}_{self.label}.pdf"))
