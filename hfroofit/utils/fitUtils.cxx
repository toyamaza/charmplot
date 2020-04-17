namespace Charm {

    void fitAndPlot(RooStats::HistFactory::Measurement& Meas, RooWorkspace* ws) {
        TFile* outFile = new TFile("out.root", "RECREATE");
        FILE* tableFile = fopen("tableFile.txt", "a");
        RooStats::HistFactory::FitModelAndPlot(Meas.GetName(), Meas.GetOutputFilePrefix(), ws, "combined", "obsData", outFile, tableFile);
    }

}
