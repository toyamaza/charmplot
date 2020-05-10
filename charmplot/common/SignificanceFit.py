from charmplot.common import canvas
from charmplot.common import utils
from charmplot.control import channel
from charmplot.control import variable
import logging
import os
import ROOT

# logging
logger = logging.getLogger(__name__)


class SigFit(object):
    
    def __init__(name: str, histogram: ROOT.TH1, label: str, output: str = ""):
        
        self.name = name
        self.variable = "Dmeson_m"
        self.h = histogram.Clone(f"{histogram.GetName()}_mass_fit")
        self.label = label
        self.output = output

        self.range = [1.75, 2.1]
        self.bkg_function = "poly"
        self.poly_order = 2
        self.SumW2Error = True

    def fitSB(self):
        x = ROOT.RooRealVar("x", "m(D^{+})", self.range[0], self.range[1], "GeV")
        x_arg_list = ROOT.RooArgList(x)
        dh = ROOT.RooDataHist("dh", "dh", x_arg_list, ROOT.RooFit.Import(self.h))

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

        ROOT.Charm.chi2Fit(bkg, dh, self.SumW2Error)

        frame = x.frame()

        # data
        dh.plotOn(frame, ROOT.RooFit.MarkerColor(ROOT.kBlack), ROOT.RooFit.Name("data"), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))

        # BKG sideband fit
        bkg.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("SB Fit"))

        # Chi2
        floating = bkg.getParameters(dh).selectByAttrib("Constant", ROOT.kFALSE).getSize()
        ndf = dh.numEntries() - floating
        chi2 = frame.chiSquare("bkgSB", "data", floating)
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