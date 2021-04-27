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

from ROOT import gROOT
atlasrootstyle_path = '/global/u1/b/bowjun/atlasrootstyle'
gROOT.SetMacroPath(os.pathsep.join([gROOT.GetMacroPath(), atlasrootstyle_path ]))
gROOT.LoadMacro("AtlasLabels.C")
gROOT.LoadMacro("AtlasStyle.C")
gROOT.LoadMacro("AtlasUtils.C")

from ROOT import (
    gRandom,
    SetAtlasStyle,
    TH1F,
)

SetAtlasStyle()

UB = 1.92999 # for D+ only
LB = 1.80001# for D+ only
SF_init = 0.3944531837982832

#function of the side bands
def SB_function(x, p):
    if x[0] < UB and x[0] > LB:
        ROOT.TF1.RejectPoint()
        return 0
    else:
        return p[0] + x[0]*p[1] + p[2]*x[0]*x[0]


if __name__ == "__main__":

    # picking which cut to look at
    var_in = "0"


    data = ["Data"]
    D_file = ROOT.TFile.Open("/global/homes/b/bowjun/minbias_charmpp/scans/multi_DpT/Lxy8_dR0.9/data/MinBias.data.root") # CHANGE TO REAL FILE PATH

    mc = ["Monte Carlo"]
    MC_file = ROOT.TFile.Open("/global/homes/b/bowjun/minbias_charmpp/scans/multi_DpT/Lxy8_dR0.9/mc/DminusKpipi.MinBias.9100011.default.root")


    data_sum = ROOT.TH1D()
    data_sum.Sumw2()
    data_sum = D_file.Get("inclusive_Dplus_per_D__Dmeson_m" + var_in).Clone("data_sum")
    #data_sum = D_file.Get("inclusive_Dplus_per_D__Dmeson_m").Clone("data_sum")
    #rebin 2 for 20 MeV bins
    data_sum.Rebin(2)

    mc_sum = ROOT.TH1D()
    mc_sum.Sumw2()
    mc_sum = MC_file.Get("inclusive_Dplus_per_D__Dmeson_m" + var_in).Clone("mc_sum")
    #rebin 2 for 20 MeV bins
    mc_sum.Rebin(2)

    # line between first and last bin = slope*(x[0] - 1.7) + (data_sum.GetBinContent(1))

    slope_data = ((data_sum.GetBinContent(data_sum.GetNbinsX()))-(data_sum.GetBinContent(1)))/((data_sum.GetXaxis().GetBinUpEdge(data_sum.GetNbinsX()))-(data_sum.GetXaxis().GetBinLowEdge(1)))
    p0_guess_data = (slope_data*(-(data_sum.GetXaxis().GetBinLowEdge(1)))) + (data_sum.GetBinContent(1))
    p1_guess_data = slope_data

    # fit to region between UB and LB using SB function
    data_bkg_fn_fit = ROOT.TF1("data_bkg_fn_fit", SB_function, 1.7, 2.2, 3)
    data_bkg_fn_fit.SetParameter(0, p0_guess_data)
    data_bkg_fn_fit.SetParameter(1, p1_guess_data)
    data_bkg_fn_fit.SetParameter(2, 0)
    data_sum.Fit(data_bkg_fn_fit, "L")
    data_bkg_fn = ROOT.TF1("data_bkg_fn", "pol2", 1.7, 2.2)
    data_bkg_fn.SetParameters(data_bkg_fn_fit.GetParameters())


    # integrate 2nd order polynomial under the curve as B
    B_data = data_bkg_fn.Integral(LB, UB) / data_sum.GetBinWidth(1)
    sqrtB_data = (B_data)**0.5
    # B_data = data_bkg_fn.Integral(1.7, 2.2)/data_bkg_fn.Integral(LB,UB)
    B_data_error = ROOT.Double()
    B_data_error = data_bkg_fn.IntegralError(LB, UB) / data_sum.GetBinWidth(1)
    # B_prob = data_bkg_fn.GetProb()
    # B_chisq = data_bkg_fn.GetChisquare()
    print ("Data Background: " + str(B_data))
    print ("Data Background '.IntegralError' error: " + str(B_data_error))
    print ("Data Background sqrtN error: " + str(sqrtB_data))
    # print ("Background fit probability: " + str(B_prob))
    # print ("Background fit chisq: " + str(B_chisq))

    T_data = data_sum.Integral(data_sum.GetXaxis().FindBin(LB), data_sum.GetXaxis().FindBin(UB))
    print ("Data Signal + Background: " + str(T_data))

    S_data = (T_data - B_data)
    print ("Data Signal: " + str(S_data))

    print()


    slope_mc = ((mc_sum.GetBinContent(mc_sum.GetNbinsX()))-(mc_sum.GetBinContent(1)))/((mc_sum.GetXaxis().GetBinUpEdge(mc_sum.GetNbinsX()))-(mc_sum.GetXaxis().GetBinLowEdge(1)))
    p0_guess_mc = (slope_mc*(-(mc_sum.GetXaxis().GetBinLowEdge(1)))) + (mc_sum.GetBinContent(1))
    p1_guess_mc = slope_mc

    mc_bkg_fn_fit = ROOT.TF1("mc_fn_fit", SB_function, 1.7, 2.2, 3)
    mc_bkg_fn_fit.SetParameter(0, p0_guess_mc)
    mc_bkg_fn_fit.SetParameter(1, p1_guess_mc)
    mc_bkg_fn_fit.SetParameter(2, 0)
    mc_sum.Fit(mc_bkg_fn_fit, "L")
    mc_bkg_fn = ROOT.TF1("mc_bkg_fn", "pol2", 1.7, 2.2)
    mc_bkg_fn.SetParameters(mc_bkg_fn_fit.GetParameters())    

    B_mc = mc_bkg_fn.Integral(LB, UB) / mc_sum.GetBinWidth(1)
    # B_mc_error = ROOT.Double()
    # B_mc_error = mc_bkg_fn.IntegralError(LB, UB) / mc_sum.GetBinWidth(1)
    # B_prob = data_bkg_fn.GetProb()
    # B_chisq = data_bkg_fn.GetChisquare()
    print ("MC Background: " + str(B_mc))
    # print ("MC Background fit error: " + str(B_mc_error))
    # print ("Background fit probability: " + str(B_prob))
    # print ("Background fit chisq: " + str(B_chisq)) 

    T_mc = mc_sum.Integral(mc_sum.GetXaxis().FindBin(LB), mc_sum.GetXaxis().FindBin(UB))
    print ("MC Signal + Background: " + str(T_mc))

    S_mc = (T_mc - B_mc)
    print ("MC Signal: " + str(S_mc))

    # data signal error
    # IntSig_mc = ROOT.Double()
    # mc_signal_error = ROOT.Double()
    # IntSig_mc = mc_sum.IntegralAndError(mc_sum.GetXaxis().FindBin(LB), mc_sum.GetXaxis().FindBin(UB), mc_signal_error)
    # mc_signal_error_Z = mc_signal_error * SF_init
    sqrtS_mc = (SF_init * S_mc)**0.5
    print("MC Signal sqrtN Error (scaled to data): " + str(sqrtS_mc))
    #print ("MC Signal '.IntegralAndError' error (scaled to data): " + str(???))
    print()

    # MC scale factor
    SF = (S_data / S_mc)
    print("Scale factor: " + str(SF))




    Z = (S_mc * SF_init) / ((sqrtB_data * sqrtB_data) + (sqrtS_mc * sqrtS_mc))**0.5
    print ("\nZ = " + str(Z))
    print ()
    



    # plot histogram with fit to background
    mc_canv = ROOT.TCanvas("mc_canv", "plot", 800, 600)
    mc_sum.GetXaxis().SetRangeUser(1.71, 2.2)
    mc_sum.GetYaxis().SetRangeUser(0, 30000)
    mc_sum.SetXTitle("MC: D+ Mass [GeV]")
    mc_sum.SetYTitle("Entries / (0.02 GeV)")
    mc_sum.SetLineColor(2)
    mc_sum.Draw("hist e")
    mc_bkg_fn.SetLineColor(2)
    mc_bkg_fn.Draw("same")
    mc_canv.Draw()
    
    # To generate the output plots:
    mc_canv.Print("mc_Lxy8_dR0.9_DpT3.0_final.png")



    d_canv = ROOT.TCanvas("d_canv", "plot", 800, 600)
    data_sum.GetXaxis().SetRangeUser(1.71, 2.2)
    data_sum.GetYaxis().SetRangeUser(0, 10000)
    data_sum.SetXTitle("D+ Mass [GeV]")
    data_sum.SetYTitle("Entries / (0.02 GeV)")
    # data_sum.SetTitle("Bkg: " + str(B) + " entries, Sig: " + str(S) + " entries")
    data_sum.SetMarkerSize(0)
    # l1 = ROOT.TLatex()
    # l1.SetTextFont(72)
    # l1.SetTextSize(28)
    # l1.DrawLatex(0.18, .85, "ATLAS Internal")
    # data_sum.SetStats(0)
    data_sum.SetLineColor(1)
    data_sum.Draw("hist e")
    data_bkg_fn.SetLineColor(2)
    data_bkg_fn.Draw("same")
    # data_bkg_fn_fit.SetLineColor(3)
    # data_bkg_fn_fit.Draw("same")
    d_canv.Draw()
    ROOT.ATLASLabel(0.20, 0.88, "Work In Progress", 1)
    ROOT.myText(0.20, 0.82, 1, "#sqrt{s} = 13 TeV, 1.637 nb^{-1}")
    ROOT.myText(0.20, 0.76, 1, "D#rightarrowK#pi#pi")
    
    # To generate the output plots:
    d_canv.Print("d_Lxy8_dR0.9_DpT3.0_final.png")