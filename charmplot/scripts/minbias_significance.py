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
SF_init = 0.36487842701502454

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


    data = ["Data"] # what does this do again..?
    D_file = ROOT.TFile.Open("/global/homes/b/bowjun/minbias_charmpp/Data3/MinBias.data.root") # CHANGE TO REAL FILE PATH

    mc = ["Monte Carlo"]
    MC_file = ROOT.TFile.Open("/global/homes/b/bowjun/minbias_charmpp/MCtests/Jan26/DminusKpipi.MinBias.9100011.default.root")


    data_sum = ROOT.TH1D()
    data_sum.Sumw2()
    data_sum = D_file.Get("inclusive_Dplus_per_D__Dmeson_m").Clone("data_sum")
    data_sum.Rebin(2)

    mc_sum = ROOT.TH1D()
    mc_sum.Sumw2()
    mc_sum = MC_file.Get("inclusive_Dplus_per_D__Dmeson_m").Clone("mc_sum")
    mc_sum.Rebin(2)



    # fit to region between UB and LB using SB function
    # Call fitting function to get integral and uncertainty
    bkg_fn_fit = ROOT.TF1("bkg_fn_fit", SB_function, 1.7, 2.2, 3)
    # bkg_fn_fit.SetParLimits(0, 50, 100)
    data_sum.Fit(bkg_fn_fit, "L")
    data_bkg_fn = ROOT.TF1("data_bkg_fn", "pol2", 1.7, 2.2)
    data_bkg_fn.SetParameters(bkg_fn_fit.GetParameters())

    # integrate 2nd order polynomial under the curve as B
    B_data = data_bkg_fn.Integral(LB, UB) / data_sum.GetBinWidth(1)
    B_data_error = ROOT.Double()
    B_data_error = data_bkg_fn.IntegralError(LB, UB) / data_sum.GetBinWidth(1)
    # B_prob = data_bkg_fn.GetProb()
    # B_chisq = data_bkg_fn.GetChisquare()
    print ("Data Background: " + str(B_data))
    print ("Data Background fit error: " + str(B_data_error))
    # print ("Background fit probability: " + str(B_prob))
    # print ("Background fit chisq: " + str(B_chisq))

    T_data = data_sum.Integral(data_sum.GetXaxis().FindBin(LB) + 1, data_sum.GetXaxis().FindBin(UB) - 1)
    print ("Data Signal + Background: " + str(T_data))

    S_data = (T_data - B_data)
    print ("Data Signal: " + str(S_data))

    print()




    mc_bkg_fn_fit = ROOT.TF1("mc_fn_fit", SB_function, 1.7, 2.0, 3)
    # mc_bkg_fn_fit.SetParLimits(???)
    mc_sum.Fit(mc_bkg_fn_fit, "L")
    mc_bkg_fn = ROOT.TF1("mc_bkg_fn", "pol2", 1.7, 2.0)
    mc_bkg_fn.SetParameters(mc_bkg_fn_fit.GetParameters())    

    B_mc = mc_bkg_fn.Integral(LB, UB) / data_sum.GetBinWidth(1)
    B_mc_error = ROOT.Double()
    B_mc_error = mc_bkg_fn.IntegralError(LB, UB) / mc_sum.GetBinWidth(1)
    # B_prob = data_bkg_fn.GetProb()
    # B_chisq = data_bkg_fn.GetChisquare()
    print ("MC Background: " + str(B_mc))
    print ("MC Background fit error: " + str(B_mc_error))
    # print ("Background fit probability: " + str(B_prob))
    # print ("Background fit chisq: " + str(B_chisq)) 

    T_mc = mc_sum.Integral(mc_sum.GetXaxis().FindBin(LB) + 1, mc_sum.GetXaxis().FindBin(UB) -1)
    print ("MC Signal + Background: " + str(T_mc))

    S_mc = (T_mc - B_mc)
    print ("MC Signal: " + str(S_mc))

    # data signal error
    IntSig_mc = ROOT.Double()
    mc_signal_error = ROOT.Double()
    IntSig_mc = mc_sum.IntegralAndError(data_sum.GetXaxis().FindBin(LB) + 1, data_sum.GetXaxis().FindBin(UB) - 1, mc_signal_error)
    mc_signal_error_Z = mc_signal_error * SF_init
    print("MC Signal Error: " + str(mc_signal_error_Z))
    print()

    # MC scale factor
    SF = (S_data / S_mc)
    print("Scale factor: " + str(SF))




    Z = S_mc * SF_init / (B_data_error * B_data_error + mc_signal_error_Z * mc_signal_error_Z)**0.5
    print ("\nZ = " + str(Z))
    print ()
    


    # plot histogram with fit to background
    mc_canv = ROOT.TCanvas("mc_canv", "plot", 800, 600)
    mc_sum.GetXaxis().SetRangeUser(1.71, 2.0)
    mc_sum.GetYaxis().SetRangeUser(0, 2500)
    mc_sum.SetXTitle("MC: D+ Mass [GeV]")
    mc_sum.SetYTitle("Entries / (0.02 GeV)")
    mc_sum.SetLineColor(2)
    mc_sum.Draw("hist e")
    mc_bkg_fn.SetLineColor(2)
    mc_bkg_fn.Draw("same")
    mc_canv.Draw()
    mc_canv.Print("mc_test.pdf")


    d_canv = ROOT.TCanvas("d_canv", "plot", 800, 600)
    data_sum.GetXaxis().SetRangeUser(1.71, 2.2)
    data_sum.GetYaxis().SetRangeUser(0, 1200)
    data_sum.SetXTitle("Data: D+ Mass [GeV]")
    data_sum.SetYTitle("Entries / (0.02 GeV)")
    # data_sum.SetTitle("Bkg: " + str(B) + " entries, Sig: " + str(S) + " entries")
    # l1 = ROOT.TLatex()
    # l1.SetTextFont(72)
    # l1.SetTextSize(28)
    # l1.DrawLatex(0.18, .85, "ATLAS Internal")
    # data_sum.SetStats(0)
    data_sum.SetLineColor(4)
    data_sum.Draw("hist e")
    data_bkg_fn.SetLineColor(4)
    data_bkg_fn.Draw("same")
    # bkg_fn_fit.SetLineColor(3)
    # bkg_fn_fit.Draw("same")
    d_canv.Draw()
    d_canv.Print("d_test.pdf")   





