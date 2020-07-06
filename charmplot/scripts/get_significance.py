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

#Global Variables
UB = 1.93
LB = 1.80
data_mc_SF = 0.9975 #0.8261
dataNminus1 = 108102.73721837193
mcNminus1 = 111457.96897521973

def SB_function(x, p):
    if x[0] < UB and x[0] > LB:
        ROOT.TF1.RejectPoint()
        return 0
    else:
        return p[0] + p[1] * x[0] + p[2] * x[0] * x[0]


if __name__ == "__main__":

    var_in = "9"

    mc_samples =   ["Wjets_emu_Rest", "Wjets_emu_Matched","Wjets_emu_NoMatch", "Top_Matched", "Top_Rest" , "Zjets_emu", "Other"]
#    mc_samples =   ["Wjets_emu_HardMisMatched", "Wjets_emu_Matched","Wjets_emu_NoMatch","Wjets_emu_Charm", "Top_Matched", "Top_Charm", "Top_HardMisMatched" , "Zjets_emu", "Other"]
    
    data_samples = ["Data"]
    file = ROOT.TFile("histograms.root")

    charges = ["OS", "SS", "OS-SS"]
#    charges = ["OS-SS"]
    leps = ["el", "mu"]
    sig_sums = []
    sig_sums_data = []
    mc_sums   = []
    data_sums = []
    
    first_hist_mc = True
    first_hist_data = True
    first_hist_sig = True
    first_hist_sig_data = True

    for c in charges:
       
        #Getting background hists
        mc_sum = ROOT.TH1F()
        mc_sum.Sumw2()
        data_sum = ROOT.TH1F()
        data_sum.Sumw2()

        for l in leps:
            # Accumulating background Histograms
            for s in mc_samples:
                if first_hist_mc:
                    mc_sum.Clear()
                    mc_sum = file.Get(s + "_2018_" + l + "_wplusd_" + c + "_Dmeson_m" + var_in).Clone()
                    first_hist_mc = False
                else:
                    hist = ROOT.TH1F()
                    print(s + "_2018_" + l + "_wplusd_" + c + "_Dmeson_m" + var_in)
                    hist = file.Get(s + "_2018_" + l + "_wplusd_" + c + "_Dmeson_m" + var_in).Clone()
                    print(s + "_2018_" + l + "_wplusd_" + c + "_Dmeson_m" + var_in)
                    mc_sum.Add(hist)
            print()    
            for s in data_samples:
                
                if first_hist_data:
                    data_sum.Clear()
                    data_sum = file.Get(s + "_2018_" + l + "_wplusd_" + c + "_Dmeson_m" + var_in).Clone()
                    #print(s + "_2018_" + l + "_wplusd_" + c + "_Dmeson_m")
                    first_hist_data = False
                else:
                    hist = file.Get(s + "_2018_" + l + "_wplusd_" + c + "_Dmeson_m" + var_in).Clone()
                    #print(s + "_2018_" + l + "_wplusd_" + c + "_Dmeson_m")
                    #print (hist)
                    data_sum.Add(hist)
            print()
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
            print("lep: "+l)
            if first_hist_sig:
                sig_sum.Clear()
                sig_sum = file.Get("Wjets_emu_Matched_2018_" + l + "_wplusd_" + c + "_Dmeson_m" + var_in).Clone()
                first_hist_sig = False
                print("el MC: " + str(sig_sum.Integral()))
            else:
                hist = ROOT.TH1F()
                hist = file.Get("Wjets_emu_Matched_2018_" + l + "_wplusd_" + c + "_Dmeson_m" + var_in).Clone()
                sig_sum.Add(hist)
                print("mu MC: " + str(hist.Integral()))
            if first_hist_sig_data:
                sig_sum_data.Clear()
                sig_sum_data = file.Get("Data_2018_" + l + "_wplusd_" + c + "_Dmeson_m" + var_in).Clone()
                first_hist_sig_data = False
                print("el MC: " + str(sig_sum_data.Integral()))
            else:
                hist_data = ROOT.TH1F()
                hist_data = file.Get("Data_2018_" + l + "_wplusd_" + c + "_Dmeson_m" + var_in).Clone()
                sig_sum_data.Add(hist_data)
                print("mu MC: " + str(hist_data.Integral()))
        first_hist_sig = True
        first_hist_sig_data = True
        sig_sums += [sig_sum]
        sig_sums_data += [sig_sum_data]
        
    print("lep sum: " + str(sig_sums[2].Integral()))
    print("lep sum: " + str(sig_sums_data[2].Integral()))


    bkg_fns = []
    for c, d, mc, sig, sig_data in zip(charges, data_sums, mc_sums, sig_sums, sig_sums_data):
        dError = ROOT.Double()
        mcError = ROOT.Double()       
        dInt = d.IntegralAndError(1,100,dError)
        mcInt = mc.IntegralAndError(1,100,mcError)
        #print ("Integral in SB for " + c +" mc:   " + str(mcInt) + " +/- " + str(mcError))

        # Call fitting function to get integral and uncertainty
        bkg_fn_fit = ROOT.TF1("bkg_fn_fit", SB_function, sig_data.GetBinLowEdge(1), sig_data.GetBinLowEdge(sig_data.GetNbinsX()) + sig_data.GetBinWidth(1), 3)
        d.Fit(bkg_fn_fit, "0")
        
        bkg_fn = ROOT.TF1("bkg_fn", "pol2", sig_data.GetBinLowEdge(1), sig_data.GetBinLowEdge(sig_data.GetNbinsX()) + sig_data.GetBinWidth(1))

        bkg_fn.SetParameters(bkg_fn_fit.GetParameters())
        bkg_fns += [bkg_fn]
        dIntSR = ROOT.Double()
        dErrorSR = ROOT.Double()
        # Events and Uncertainty in SB region 
        dIntSR = bkg_fn.Integral(LB, UB) / d.GetBinWidth(1)
        dErrorSR = bkg_fn.IntegralError(LB, UB) / d.GetBinWidth(1)
        # Events in Matched Peak Region
        dIntSig = ROOT.Double()
        dIntSig = sig.Integral(sig.GetXaxis().FindBin(LB) + 1 , sig.GetXaxis().FindBin(UB) - 1)
        # Events in data SR used to normalize MC
        dIntSig_data = ROOT.Double()
        dIntSig_data_err = ROOT.Double()
        dIntSig_data = sig_data.IntegralAndError(sig_data.GetXaxis().FindBin(LB) + 1 , sig_data.GetXaxis().FindBin(UB) - 1, dIntSig_data_err)
        print ("just dIntSig_integral: " + str(sig_data.Integral(sig_data.GetXaxis().FindBin(LB) + 1 , sig_data.GetXaxis().FindBin(UB) - 1)))
        #Replace uncertainty on Signal Error for OS - SS case
        if "OS-SS" in c:
            #OS Signal Error
            OS_err2 = sig_sums[0].Integral(sig.GetXaxis().FindBin(LB) + 1, sig.GetXaxis().FindBin(UB) - 1)
            SS_err2 = sig_sums[1].Integral(sig.GetXaxis().FindBin(LB) + 1, sig.GetXaxis().FindBin(UB) - 1)
            dIntSig_data_err = (OS_err2 + SS_err2)**.5 * data_mc_SF
        # Significance where Z = N events in signal peak / uncertainty on bkg sideband
        
        Z = dIntSig * data_mc_SF / (dErrorSR * dErrorSR + dIntSig_data_err * dIntSig_data_err)**0.5

        canv = ROOT.TCanvas(c,c,800,600)
        d.GetXaxis().SetRangeUser(1.71,2.2)
        d.GetYaxis().SetRangeUser(d.GetMinimum(),d.GetMaximum()*1.3)
        d.GetXaxis().SetTitle("D Mass [GeV]")
        l1 = ROOT.TLatex()
        l1.SetTextFont(72)
        l1.SetTextSize(28)
        l1.DrawLatex(0.18, .85, "ATLAS Internal")
        d.Draw("")
        bkg_fn.SetLineColor(2)
        bkg_fn.Draw("same")
        canv.Print(c+"_test.pdf")

        bkgHist = bkg_fn.CreateHistogram()
        print("Nbins bkg, sig: " + str(bkgHist.GetNbinsX()) + ", " + str(sig_data.GetNbinsX()))
        print("0 bin bkg, sig: " + str(bkgHist.GetBinCenter(100)) + ", " + str(sig_data.GetBinCenter(100)))
        sigMinusBkg = sig_data
        sigMinusBkg.Add(bkgHist,-1)
        sigMinusBkg.GetXaxis().SetRangeUser(1.76,2.0)
        sigMinusBkg.GetYaxis().SetRangeUser(sigMinusBkg.GetMinimum(), sigMinusBkg.GetMaximum()*1.3)
        sigMinusBkg.GetXaxis().SetTitle("D Mass [GeV]")
        sigMinusBkg.Draw()
        sig.Scale(data_mc_SF)
        sig.SetLineColor(2)
        sig.Draw("histsame")
        canv.Print(c+"_test2.pdf")

        print( "Chi2/NDF:  " + str(bkg_fn_fit.GetChisquare()/bkg_fn_fit.GetNDF()))
        #print ("Integral in Peak for Signal + Background " + c + " data: " + str(dInt) + " +/- " + str(dError))
        #print ("Integral in SR for fitted bkg to " + c + " data:         " + str(dIntSR) + " +/- " + str(dErrorSR))
        print ("Integral in SR for Scaled Matched MC in " + c  + "              " + str(dIntSig * data_mc_SF))
        #print ("Integral for data - bkg fn in peak: " + str(sigMinusBkg.Integral(sig.GetXaxis().FindBin(LB) + 1, sig.GetXaxis().FindBin(UB) - 1)))
        print ("dIntSig_data - dIntSR: " + str(dIntSig_data - dIntSR))
        print ("Significance for " + c + "                               " + str(Z))
        print ("Data Eff: " + str((dIntSig_data - dIntSR)/dataNminus1))
        print ("MC   Eff: " + str((dIntSig * data_mc_SF)/mcNminus1))
        print ("Normalization for data/MC in peak region vs Baseline SF :               " + str((dIntSig_data - dIntSR)/(dIntSig)) + " : " + str(data_mc_SF))