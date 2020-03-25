namespace Charm {

    RooFitResult* chi2Fit(RooAbsPdf* pdf, RooDataHist* dh, bool SumW2Error) {
        return pdf->chi2FitTo(*dh, RooFit::Minimizer("Minuit2", "Migrad"),
                              RooFit::SumW2Error(SumW2Error),
                              RooFit::Save(),
                              RooFit::PrintLevel(9));
    }

    RooFitResult* likelihoodFit(RooAbsPdf* pdf, RooDataHist* dh, bool SumW2Error) {
        return pdf->fitTo(*dh, RooFit::Minimizer("Minuit2", "Migrad"),
                          RooFit::SumW2Error(SumW2Error),
                          RooFit::Save(),
                          RooFit::PrintLevel(9));
    }

}
