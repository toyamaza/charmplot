import os
import ROOT
# import numpy as np
# import scipy.integrate as integrate
# from scipy.integrate import quad
from charmplot.common import utils

# ATLAS Style
# dirname = os.path.join(os.path.dirname(__file__), "../../atlasrootstyle")
# ROOT.gROOT.SetBatch(True)
# ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
# ROOT.SetAtlasStyle()

UB = 1.93 # for D+ only
LB = 1.80# for D+ only

#function of the side bands
def SB_function(x, p):
    if x[0] < UB and x[0] > LB:
        ROOT.TF1.RejectPoint()
        return 0
    else:
        return p[0] + x[0]*p[1] + p[2]*x[0]*x[0]

# ROOT.gInterpreter.ProcessLine("""
# Double_t SB_function(Double_t *x, Double_t *par)
# {
#     if (x[0] > 1.80 && x[0] < 1.93) {
#       TF1::RejectPoint();
#       return 0;
#    }
#    return par[0] + par[1]*x[0] + par[2]*x[0]*x[0];
# }
# """)

if __name__ == "__main__":

    # picking which cut to look at
    #var_in = ""

    # loading in data
    data = ["Data"] # what does this do again..?
    file = ROOT.TFile.Open("/global/u1/b/bowjun/minbias_charmpp//Data2/Dplus_all/MinBias.data.00267599.root") # CHANGE TO REAL FILE PATH

    # read in a root histogram file of var_in
    # data_sum = []
    
    # data_sum = ROOT.TH1F()
    # data_sum.Sumw2()
    
    # data_sum.Clear()
    data_sum = file.Get("inclusive_Dplus_per_D/inclusive_Dplus_per_D__Dmeson_m").Clone()

    data_sum.Rebin(0)

    # fit to region between UB and LB using SB function
    # Call fitting function to get integral and uncertainty
    bkg_fn_fit = ROOT.TF1("bkg_fn_fit", SB_function, 1.7, 2.2, 3)
    # bkg_fn_fit.SetParLimits(0, 50, 100)
    data_sum.Fit(bkg_fn_fit, "0")

    bkg_fn = ROOT.TF1("bkg_fn", "pol2", 1.7, 2.2)
    bkg_fn.SetParameters(bkg_fn_fit.GetParameters())

    # integrate 2nd order polynomial under the curve as B
    B = bkg_fn.Integral(LB, UB) / data_sum.GetBinWidth(1)
    B_error = bkg_fn.IntegralError(LB, UB) / data_sum.GetBinWidth(1)
    # B_prob = bkg_fn.GetProb()
    # B_chisq = bkg_fn.GetChisquare()
    print ("Background: " + str(B))
    print ("Background fit error: " + str(B_error))
    # print ("Background fit probability: " + str(B_prob))
    # print ("Background fit chisq: " + str(B_chisq))

    # integrate full peak with background under as T total signal + background
    T = data_sum.Integral(data_sum.GetXaxis().FindBin(LB) + 1, data_sum.GetXaxis().FindBin(UB) - 1)
    print ("Signal + Background: " + str(T))

    # subtract B from T to get S (the numerator of Z)
    S = (T - B)
    print ("Signal: " + str(S))

    # NEXT STEP? find unc in background, and get Z = S / sigmaB


    # plot histogram with fit to background

    canv = ROOT.TCanvas("canvas", "plot", 800, 600)
    data_sum.GetXaxis().SetRangeUser(1.71, 2.2)
    data_sum.GetYaxis().SetRangeUser(0, 720)
    data_sum.SetXTitle("D+ Mass [GeV]")
    data_sum.SetYTitle("Entries / (0.020 GeV)")
    data_sum.SetTitle("Bkg: " + str(B) + " entries, Sig: " + str(S) + " entries")
    # l1 = ROOT.TLatex()
    # l1.SetTextFont(72)
    # l1.SetTextSize(28)
    # l1.DrawLatex(0.18, .85, "ATLAS Internal")
    # data_sum.SetStats(0)
    data_sum.SetLineColor(1)
    data_sum.Draw("hist e")
    bkg_fn.SetLineColor(2)
    bkg_fn.Draw("same")
    # bkg_fn_fit.SetLineColor(3)
    # bkg_fn_fit.Draw("same")

    canv.Print("data_fit_00267599.pdf")   

