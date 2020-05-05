#!/usr/bin/env python
import os
import ROOT
from charmplot.common import utils

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.SetAtlasStyle()

#Global Variables
UB = 1.95
LB = 1.79

def SB_function(x, p):
    if x[0] < UB and x[0] > LB:
        ROOT.TF1.RejectPoint()
        return 0
    else:
        return p[0] + p[1] * x[0] + p[2] * x[0] * x[0]


if __name__ == "__main__":

    mc_samples =   ["Wjets_emu_Rest", "Wjets_emu_Matched", "Wjets_tau", "Top", "Diboson", "Zjets"]
    data_samples = ["Data"]
    file = ROOT.TFile("histograms.root")

    charges = ["OS", "SS", "OS-SS"]
    leps = ["el", "mu"]
    sig_sums = []
    sig_sums_data = []
    mc_sums   = []
    data_sums = []
    
    first_hist_mc = True
    first_hist_data = True
    first_hist_sig = True

    for c in charges:
        # Get Matched MC hist for signal region
        sig_sum = ROOT.TH1F()
        sig_sum.Sumw2()
        sig_sum_data = ROOT.TH1F()
        sig_sum_data.Sumw2()
        #Getting background hists
        mc_sum = ROOT.TH1F()
        mc_sum.Sumw2()
        data_sum = ROOT.TH1F()
        data_sum.Sumw2()
        for l in leps:
            if first_hist_sig:
                sig_sum.Clear()
                sig_sum = file.Get("Wjets_emu_Matched_2018_" + l + "_wplusd_" + c + "_Dmeson_m")
                sig_sum_data = file.Get("Data_2018_" + l + "_wplusd_" + c + "_Dmeson_m")
                first_hist_sig = False
            else:
                hist = ROOT.TH1F()
                hist_data = ROOT.TH1F()
                hist = file.Get("Wjets_emu_Matched_2018_" + l + "_wplusd_" + c + "_Dmeson_m")
                hist_data = file.Get("Data_2018_" + l + "_wplusd_" + c + "_Dmeson_m")
                sig_sum.Add(hist)
                sig_sum_data.Add(hist_data)
            # Accumulating background Histograms
            for s in mc_samples:
                if first_hist_mc:
                    mc_sum.Clear()
                    mc_sum = file.Get(s + "_2018_" + l + "_wplusd_" + c + "_Dmeson_m")
                    first_hist_mc = False
                else:
                    hist = ROOT.TH1F()
                    hist = file.Get(s + "_2018_" + l + "_wplusd_" + c + "_Dmeson_m")
                    mc_sum.Add(hist)
            print()    
            for s in data_samples:
                if first_hist_data:
                    data_sum = file.Get(s + "_2018_" + l + "_wplusd_" + c + "_Dmeson_m")
                    first_hist_data = False
                else:
                    hist = file.Get(s + "_2018_" + l + "_wplusd_" + c + "_Dmeson_m")
                    #print (hist)
                    data_sum.Add(hist)
            print()
        first_hist_mc = True
        first_hist_data = True
        first_hist_sig = True
        sig_sums += [sig_sum]
        sig_sums_data += [sig_sum_data]
        mc_sums += [mc_sum]
        data_sums += [data_sum]
        print(mc_sums)
        print(data_sums)


    bkg_fns = []
    for c, d, mc, sig, sig_data in zip(charges, data_sums, mc_sums, sig_sums, sig_sums_data):
        dError = ROOT.Double()
        mcError = ROOT.Double()       
        dInt = d.IntegralAndError(1,100,dError)
        mcInt = mc.IntegralAndError(1,100,mcError)
        #print ("Integral in SB for " + c +" mc:   " + str(mcInt) + " +/- " + str(mcError))

        # Call fitting function to get integral and uncertainty
        bkg_fn_fit = ROOT.TF1("bkg_fn_fit", SB_function, 1.7, 2.2, 3)
        d.Fit(bkg_fn_fit, "0")
        
        bkg_fn = ROOT.TF1("bkg_fn", "pol2", d.GetBinLowEdge(1), d.GetBinLowEdge(d.GetNbinsX()) + d.GetBinWidth(1))

        bkg_fn.SetParameters(bkg_fn_fit.GetParameters())
        bkg_fns += [bkg_fn]
        dIntSR = ROOT.Double()
        dErrorSR = ROOT.Double()
        # Events and Uncertainty in SB region 
        dIntSR = bkg_fn.Integral(1.79, 1.95) / d.GetBinWidth(1)
        dErrorSR = bkg_fn.IntegralError(1.79, 1.95) / d.GetBinWidth(1)
        # Events in Matched Peak Region
        dIntSig = ROOT.Double()
        dIntSig = sig.Integral(sig.GetXaxis().FindBin(LB) + 1 , sig.GetXaxis().FindBin(UB) - 1)
        # Events in data SR used to normalize MC
        dIntSig_data = ROOT.Double()
        dIntSig_data_err = ROOT.Double()
        dIntSig_data = sig_data.IntegralAndError(sig_data.GetXaxis().FindBin(LB) + 1 , sig.GetXaxis().FindBin(UB) - 1, dIntSig_data_err)
        # Significance where Z = N events in signal peak / uncertainty on bkg sideband
        Z = dIntSig / (dErrorSR * dErrorSR + dIntSig_data_err * dIntSig_data_err)**0.5

        canv = ROOT.TCanvas(c,c,800,600)
        d.GetXaxis().SetRangeUser(1.7,2.2)
        d.Draw("")
        bkg_fn.SetLineColor(2)
        bkg_fn.Draw("same")
        canv.Print(c+"_test.pdf")

        bkgHist = bkg_fn.GetHistogram()
        sigMinusBkg = sig_data
        sigMinusBkg.Add(bkgHist,-1)
        sigMinusBkg.GetXaxis().SetRangeUser(1.71,2.2)
        sigMinusBkg.Draw()
        sig.Scale((dIntSig_data - dIntSR)/(dIntSig))
        sig.Draw("histsame")
        canv.Print(c+"_test2.pdf")

        print( "param 0:  " + str(bkg_fn.GetParameter(0)))
        print ("Integral in Peak for Signal + Background " + c + " data: " + str(dInt) + " +/- " + str(dError))
        print ("Integral in SR for fitted bkg to " + c + " data:         " + str(dIntSR) + " +/- " + str(dErrorSR))
        print ("Integral in SR for Matched MC in " + c  + "              " + str(dIntSig))
        print ("Significance for " + c + "                               " + str(Z))
        print ("Normalization for data/MC in peak region :               " + str((dIntSig_data - dIntSR)/(dIntSig)))