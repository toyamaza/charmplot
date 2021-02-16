# flake8: noqa
#!/usr/bin/env python
import os
import ROOT
from charmplot.common import utils

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.SetAtlasStyle()

# Global Variables
UB = 149.
LB = 143.
LimLeft = 139.5
LimRight = 175.0
#run once then set these numbers from the first line of output
data_mc_SF = 1.1748586732455397  # 0.8261
NumberEvtNminus1 = 89671.05704262384

#fitting function which excludes points from between LB and UB
def SB_function(x, p):
    if x[0] < UB and x[0] > LB:
        ROOT.TF1.RejectPoint()
        return 0
    else:
        return p[0] + p[1] * ROOT.TMath.Log(x[0] - p[2])

#Drawing and histogramming function with same form as SB_function but includes the peak region
def full_bkg_function(x,p):
    if x[0] < LimLeft or x[0] > LimRight:
        ROOT.TF1.RejectPoint()
        return 0
    else:
        return p[0] + p[1] * ROOT.TMath.Log(x[0] - p[2])


if __name__ == "__main__":

    var_in = ""

#Names of MC samples in histograms.root file
    mc_samples = ["Multijet_MatrixMethod","Wjets_emu_Rest", "Wjets_emu_Matched", "Top_Rest", "Top_Matched", "Zjets_emu", "Other"]

    #Names of Data samples in histograms.root
    data_samples = ["Data"]

    #input file in same directory as where you're running
    file = ROOT.TFile("histograms.root")

#Need OS/SS separately for uncertainty calc
    charges = ["OS", "SS", "OS-SS"]
    leps = ["el", "mu"]

#loop over all 10 variations
    for i in range(0,10):
        var_in = str(i)
        sig_sums = []
        sig_sums_data = []
        mc_sums = []
        data_sums = []

        first_hist_mc = True
        first_hist_data = True
        first_hist_sig = True
        first_hist_sig_data = True

        for c in charges:
            # Getting background hists
            mc_sum = ROOT.TH1F()
            mc_sum.Sumw2()
            data_sum = ROOT.TH1F()
            data_sum.Sumw2()

            for l in leps:
                # Accumulating background Histograms
                for s in mc_samples:
                    #All MC (Used for closure tests, but not currently implemented, can leave in without too much concern)
                    if first_hist_mc:
                        mc_sum.Clear()
                        mc_sum = file.Get(s + "_" + c +"_" + l + "_SR_Dstar_Dmeson_mdiff" + var_in).Clone()
                        first_hist_mc = False
                    else:
                        hist = ROOT.TH1F()
                        hist = file.Get(s + "_" + c +"_" + l + "_SR_Dstar_Dmeson_mdiff" + var_in).Clone()
                        mc_sum.Add(hist)
                for s in data_samples:
                    if first_hist_data:
                        data_sum.Clear()
                        data_sum = file.Get(s + "_" + c +"_" + l + "_SR_Dstar_Dmeson_mdiff" + var_in).Clone()
                        first_hist_data = False
                    else:
                        hist = file.Get(s + "_" + c +"_" + l + "_SR_Dstar_Dmeson_mdiff" + var_in).Clone()
                        data_sum.Add(hist)
            first_hist_mc = True
            first_hist_data = True
            mc_sums += [mc_sum]
            data_sums += [data_sum]

        for c in charges:
            # Get Matched MC hist for signal region
            sig_sum = ROOT.TH1F()
            sig_sum.Sumw2()
            sig_sum_data = ROOT.TH1F()
            sig_sum_data.Sumw2()

            for l in leps:
                if first_hist_sig:
                    sig_sum.Clear()
                    sig_sum = file.Get("Wjets_emu_Matched_" + c + "_" + l + "_SR_Dstar_Dmeson_mdiff" + var_in).Clone()
                    first_hist_sig = False
                else:
                    hist = ROOT.TH1F()
                    hist = file.Get("Wjets_emu_Matched_" + c + "_" + l + "_SR_Dstar_Dmeson_mdiff" + var_in).Clone()
                    sig_sum.Add(hist)
                if first_hist_sig_data:
                    sig_sum_data.Clear()
                    sig_sum_data = file.Get("Data_" + c + "_" + l + "_SR_Dstar_Dmeson_mdiff" + var_in).Clone()
                    first_hist_sig_data = False
                else:
                    hist_data = ROOT.TH1F()
                    hist_data = file.Get("Data_" + c + "_" + l + "_SR_Dstar_Dmeson_mdiff" + var_in).Clone()
                    sig_sum_data.Add(hist_data)
            first_hist_sig = True
            first_hist_sig_data = True
            sig_sums += [sig_sum]
            sig_sums_data += [sig_sum_data]

        bkg_fns = []
        for c, d, sig, sig_data in zip(charges, data_sums, sig_sums, sig_sums_data):

            # Call fitting function to get integral and uncertainty
            bkg_fn_fit = ROOT.TF1("bkg_fn_fit", SB_function, 135., 180., 3)

            # Use function below if one or two fits seem stuck in a weird local min
            #bkg_fn_fit.SetParameter(0,6.45528e+02)
            #bkg_fn_fit.SetParameter(1,1.86633e+02)
            #bkg_fn_fit.SetParameter(2,1.39683e+02)
            d.Fit(bkg_fn_fit, "Q0","",LimLeft, LimRight)

            bkg_fn = ROOT.TF1("bkg_fn", full_bkg_function, 135., 180., 3)

            bkg_fn.SetParameter(0,bkg_fn_fit.GetParameter(0))
            bkg_fn.SetParameter(2,bkg_fn_fit.GetParameter(2))
            bkg_fn.SetParameter(1,bkg_fn_fit.GetParameter(1))
           
            bkg_fns += [bkg_fn]
            dFittedBkgInt = ROOT.Double()
            dErrorFittedBkgInt = ROOT.Double()
            # Events and Uncertainty in SB region
            dFittedBkgInt = bkg_fn.Integral(LB, UB) / d.GetBinWidth(1)
            dErrorFittedBkgInt = bkg_fn.IntegralError(LB, UB) / d.GetBinWidth(1)
            # Events in Matched Peak Region
            dMcSigInt = ROOT.Double()
            dMcSigInt = sig.Integral(sig.GetXaxis().FindBin(LB) + 1, sig.GetXaxis().FindBin(UB) - 1)
            # Events in data SR used to normalize MC
            dDataSigInt = ROOT.Double()
            dDataSigIntErr = ROOT.Double()
            dDataSigInt = sig_data.IntegralAndError(sig_data.GetXaxis().FindBin(LB) + 1, sig_data.GetXaxis().FindBin(UB) - 1, dDataSigIntErr)

            # Replace uncertainty on Signal Error for OS - SS case
            if "OS-SS" in c:
                # OS Signal Error
                OS_err2 = sig_sums[0].Integral(sig.GetXaxis().FindBin(LB) + 1, sig.GetXaxis().FindBin(UB) - 1)
                SS_err2 = sig_sums[1].Integral(sig.GetXaxis().FindBin(LB) + 1, sig.GetXaxis().FindBin(UB) - 1)
                dDataSigIntErr = (OS_err2 + SS_err2)**.5 * data_mc_SF

            # Significance where Z = N events in signal peak / uncertainty on bkg sideband
            Z = dMcSigInt * data_mc_SF / (dErrorFittedBkgInt * dErrorFittedBkgInt + dDataSigIntErr * dDataSigIntErr)**0.5

            canv = ROOT.TCanvas(c, c, 800, 600)
            d.GetXaxis().SetRangeUser(LimLeft, LimRight)
            d.GetYaxis().SetRangeUser(d.GetMinimum(), d.GetMaximum() * 1.3)
            d.GetXaxis().SetTitle("D* - D^0 Mass Diff [GeV]")
            l1 = ROOT.TLatex()
            l1.SetTextFont(72)
            l1.SetTextSize(28)
            l1.DrawLatex(0.18, .85, "ATLAS Internal")
            d.Draw("")
            bkg_fn.SetLineColor(2)
            bkg_fn.Draw("same")
            canv.Print(c + "_bkgFit" + var_in + ".pdf")

            #Make sure your histograms have 100 bins 
            bkgHist = bkg_fn.CreateHistogram()
            #print("Nbins bkg, sig    : " + str(bkgHist.GetNbinsX()) + ", " + str(sig_data.GetNbinsX()))
            #print("0 bin bkg, sig    : " + str(bkgHist.GetBinCenter(100)) + ", " + str(sig_data.GetBinCenter(100)))
            #print("bin width bkg, sig: " + str(bkgHist.GetBinWidth(1)) + ", " + str(sig_data.GetBinWidth(1)))
            sigMinusBkg = sig_data
            sigMinusBkg.Add(bkgHist, -1)
            sigMinusBkg.GetXaxis().SetRangeUser(LimLeft, LimRight)
            sigMinusBkg.GetYaxis().SetRangeUser(sigMinusBkg.GetMinimum(), sigMinusBkg.GetMaximum() * 1.3)
            sigMinusBkg.GetXaxis().SetTitle("D* - D^0 Mass Diff [GeV]")
            sigMinusBkg.Draw()
            sig.Scale(data_mc_SF)
            sig.SetLineColor(2)
            sig.Draw("histsame")
            canv.Print(c + "_Data_vs_MC" + var_in + ".pdf")

            if c is "OS-SS":
                #print the lines of the table so it is easily copy/pasted to excel or similar
                print(str((dDataSigInt - dFittedBkgInt) / (dMcSigInt)) + #MC scale factor (# data signal / # MC signal)
                ", " + str(dDataSigInt - dFittedBkgInt) + # Number of events in data
                ", " + str(dMcSigInt* data_mc_SF) +  # (scaled) number of events in MC
                ", " + str((dDataSigInt - dFittedBkgInt) / NumberEvtNminus1) +  # Efficiency in Data wrt N-1 events
                ", " + str((dMcSigInt * data_mc_SF) / NumberEvtNminus1) + # Efficiency in MC wrt N-1 events
                ", " + str(((dDataSigInt - dFittedBkgInt) / NumberEvtNminus1)/((dMcSigInt * data_mc_SF) / NumberEvtNminus1))+ # ratio of data Efficiency to MC efficiency
                ", " + str(Z)) #approximate signal strength
