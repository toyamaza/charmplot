# flake8: noqa
#!/usr/bin/env python
import os
import ROOT
from charmplot.common import utils
from array import array

# ATLAS Styleget_significance.py
dirname = os.path.join(os.path.dirname(__file__), "../../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.SetAtlasStyle()

# Global Variables
UB = 1.95 #Mass GeV upper bound signal region
LB = 1.79 #Mass GeV lower bound signal region

def SB_function(x, p):
    if x[0] < UB and x[0] > LB:
        ROOT.TF1.RejectPoint()
        return 0
    else:
        return p[0] + p[1] * x[0] + p[2] * x[0] * x[0]

def chi2_function(x):
    return ROOT.Math.chisquared_pdf(x[0],3.)

if __name__ == "__main__":

    var_in = ""
    isFwd = True
    xlabel = ""
    cut_vals = []
    cut_effs_err = []
    cut_effs = []
    var_test ="chi2"
    #var_test ="Lxy"
    #var_test = "ptcone40_pt"
    #var_test = "impact_sig"

    #ntags = ["_0tag", "_1tag", "_2tag"]
    ntags = ["_2tag"]

    if var_test is "Lxy":
        xlabel = "L_{xy}"
        isFwd = False
        cut_vals = [0.0, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6]
        cut_effs = [1.000, 1.000, 0.964, 0.932, 0.896, 0.865, 0.838, 0.811, 0.784, 0.760]
        cut_effs_err = []
    elif var_test is "chi2":
        xlabel = "D \chi ^{2}"
        cut_vals = [10000, 20, 18, 15, 12, 10, 8, 7, 6, 5]
        cut_effs_err = []
        cut_effs = [1.000, 0.967, 0.956, 0.951, 0.935, 0.912, 0.878, 0.846, 0.808, 0.742]
    elif var_test is "ptcone40_pt":
        cut_vals = [100, 4.5, 4.0, 3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5]
        cut_effs = [1.000, 0.988, 0.988, 0.987, 0.987, 0.985, 0.982, 0.970, 0.917, 0.709]
        cut_effs_err = []
        xlabel = "Isolation"
    elif var_test is "impact_sig":
        cut_vals = [7, 6, 5, 4, 3]
        cut_effs = [1.000, 1.000, 0.997, 0.980, 0.922]
        cut_effs_err = [0.0041, 0.0041, 0.0041, 0.0041, 0.0039]
        xlabel = "\sigma _{3D}"

    cut_vals_c = array('d',cut_vals)
    cut_effs_c = array('d',cut_effs)
    cut_graph_errors = ROOT.TGraph(len(cut_vals),cut_vals_c, cut_effs_c)


    mc_samples_matched = ["Wjets_emu_Matched","Top_Matched"]
    mc_samples = ["Wjets_emu_Rest"]#, "Wjets_emu_Matched", "Top_Matched", "Top_Rest", "Zjets_emu", "Other", "Multijet_MatrixMethod"]

    data_samples = ["Data"]
    file = ROOT.TFile("histograms.root")

    charges = ["OS", "SS", "OS-SS"]
    leps = ["el", "mu"]
    years = ["2018"]

    mc_sums = []
    mc_sums_Nm1 = []
    mc_sums_all = []
    data_sums = []
    data_sums_Nm1 = []
    mc_sums_Chi2InSB = []
    mc_sums_all_Chi2InSB = []
    mc_sums_all_Chi2OutSB_R = []
    mc_sums_all_Chi2OutSB_L = []
    data_sums_Chi2InSB = []
    data_sums_Chi2OutSB_R = []
    data_sums_Chi2OutSB_L = []
    

    first_hist_mc = True
    first_hist_mc_all = True
    first_hist_data = True

    for c in charges:

        # Getting background hists
        mc_sum = ROOT.TH1F()
        mc_sum.Sumw2()
        mc_sum_Nm1 = ROOT.TH1F()
        mc_sum_Nm1.Sumw2()
        mc_sum_all = ROOT.TH1F()
        mc_sum_all.Sumw2()
        data_sum = ROOT.TH1F()
        data_sum.Sumw2()
        data_sum_Nm1 = ROOT.TH1F()
        data_sum_Nm1.Sumw2()
        mc_sum_Chi2InSB = ROOT.TH1F()
        mc_sum_all_Chi2InSB = ROOT.TH1F()
        mc_sum_all_Chi2OutSB_R = ROOT.TH1F()
        mc_sum_all_Chi2OutSB_L = ROOT.TH1F()
        data_sum_Chi2InSB = ROOT.TH1F()
        data_sum_Chi2OutSB_R = ROOT.TH1F()
        data_sum_Chi2OutSB_L = ROOT.TH1F()
        mc_sum_Chi2InSB.Sumw2() 
        mc_sum_all_Chi2InSB.Sumw2() 
        mc_sum_all_Chi2OutSB_R.Sumw2() 
        mc_sum_all_Chi2OutSB_L.Sumw2() 
        data_sum_Chi2InSB.Sumw2() 
        data_sum_Chi2OutSB_R.Sumw2()
        data_sum_Chi2OutSB_L.Sumw2()

        for l in leps:
            for t in ntags:
                for y in years:
                    # Accumulating background Histograms
                    for s in mc_samples_matched:
                        if first_hist_mc:
                            mc_sum.Clear()
                            mc_sum = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_m" + var_in).Clone()
                            mc_sum_Nm1.Clear()
                            mc_sum_Nm1 = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_m_" + var_test  + var_in).Clone()
                            mc_sum_Chi2InSB.Clear()
                            mc_sum_Chi2InSB = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_" + var_test + "_inSB" + var_in).Clone()
                            first_hist_mc = False
                        else:
                            hist = ROOT.TH1F()
                            hist2 = ROOT.TH1F()
                            hist21 = ROOT.TH1F()
                            print(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_m" + var_in)
                            hist = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_m" + var_in).Clone()
                            hist21 = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_m_" + var_test + var_in).Clone()
                            hist2 = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_" + var_test + "_inSB" + var_in).Clone()
                            print(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_m" + var_in)
                            mc_sum.Add(hist)
                            mc_sum_Nm1.Add(hist21)
                            mc_sum_Chi2InSB.Add(hist2)
                    print()
                    for s in mc_samples:
                        if first_hist_mc_all:
                            mc_sum_all.Clear()
                            mc_sum_all = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_m" + var_in).Clone()
                            mc_sum_all_Chi2InSB.Clear()
                            mc_sum_all_Chi2InSB = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_" + var_test + "_inSB" + var_in).Clone()
                            mc_sum_all_Chi2OutSB_R.Clear()
                            mc_sum_all_Chi2OutSB_L.Clear()
                            mc_sum_all_Chi2OutSB_R = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_" + var_test + "_outSB_R" + var_in).Clone()
                            mc_sum_all_Chi2OutSB_L = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_" + var_test + "_outSB_L" + var_in).Clone()
                            first_hist_mc_all = False
                        else:
                            hist_all = ROOT.TH1F()
                            hist2_all = ROOT.TH1F()
                            hist3_all = ROOT.TH1F()
                            hist4_all = ROOT.TH1F()
                            print(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_m" + var_in)
                            hist_all = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_m" + var_in).Clone()
                            hist2_all = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_" + var_test + "_inSB" + var_in).Clone()
                            hist3_all = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_" + var_test + "_outSB_R" + var_in).Clone()
                            hist4_all = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_" + var_test + "_outSB_L" + var_in).Clone()
                            print(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_m" + var_in)
                            mc_sum_all.Add(hist_all)
                            mc_sum_all_Chi2InSB.Add(hist2_all)
                            mc_sum_all_Chi2OutSB_R.Add(hist3_all)
                            mc_sum_all_Chi2OutSB_L.Add(hist4_all)
                    print()

                    for s in data_samples:
                        if first_hist_data:
                            data_sum.Clear()
                            data_sum_Nm1.Clear()
                            data_sum_Chi2InSB.Clear()
                            data_sum_Chi2OutSB_R.Clear()
                            data_sum_Chi2OutSB_L.Clear()
                            data_sum = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_m" + var_in).Clone()
                            data_sum_Nm1 = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_m_" + var_test + var_in).Clone()
                            data_sum_Chi2InSB = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_" + var_test + "_inSB" + var_in).Clone()
                            data_sum_Chi2OutSB_R = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_" + var_test + "_outSB_R" + var_in).Clone()
                            data_sum_Chi2OutSB_L = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_" + var_test + "_outSB_L" + var_in).Clone()
                            #print(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_m")
                            first_hist_data = False
                        else:
                            hist3 = ROOT.TH1F()
                            hist4 = ROOT.TH1F()
                            hist5 = ROOT.TH1F()
                            hist = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_m" + var_in).Clone()
                            hist5 = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_m_" + var_test + var_in).Clone()
                            hist2 = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_" + var_test + "_inSB" + var_in).Clone()
                            hist3 = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_" + var_test + "_outSB_R" + var_in).Clone()
                            hist4 = file.Get(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_" + var_test + "_outSB_L" + var_in).Clone()
                            #print(s + "_" + c + "_" + y +"_" + l + t + "_SR_Dplus_Dmeson_m")
                            #print (hist)
                            data_sum.Add(hist)
                            data_sum_Nm1.Add(hist5)
                            data_sum_Chi2InSB.Add(hist2) 
                            data_sum_Chi2OutSB_R.Add(hist3)
                            data_sum_Chi2OutSB_L.Add(hist4)
            print()
        first_hist_mc = True
        first_hist_data = True
        first_hist_mc_all = True
        mc_sums += [mc_sum]
        mc_sums_Nm1 += [mc_sum_Nm1]
        mc_sums_all += [mc_sum_all]
        data_sums += [data_sum]
        data_sums_Nm1 += [data_sum_Nm1]
        mc_sums_Chi2InSB += [mc_sum_Chi2InSB]
        mc_sums_all_Chi2InSB += [mc_sum_all_Chi2InSB]
        mc_sums_all_Chi2OutSB_R += [mc_sum_all_Chi2OutSB_R]
        mc_sums_all_Chi2OutSB_L += [mc_sum_all_Chi2OutSB_L]
        data_sums_Chi2InSB += [data_sum_Chi2InSB]
        data_sums_Chi2OutSB_R += [data_sum_Chi2OutSB_R]
        data_sums_Chi2OutSB_L += [data_sum_Chi2OutSB_L]

    bkg_fns = []
    OS_err = []
    SS_err = []
    OS_mc_err = []
    SS_mc_err = []
    for c, d, d_Nm1, mc, mc_Nm1, mc_all, mc_chi2, mc_all_varTest_inSB, mc_all_varTest_outSB_R, mc_all_varTest_outSB_L, data_chi2InSB, data_chi2OutSB_R, data_chi2OutSB_L in zip(charges, data_sums, data_sums_Nm1, mc_sums, mc_sums_Nm1, mc_sums_all, mc_sums_Chi2InSB, mc_sums_all_Chi2InSB, mc_sums_all_Chi2OutSB_R, mc_sums_all_Chi2OutSB_L, data_sums_Chi2InSB, data_sums_Chi2OutSB_R, data_sums_Chi2OutSB_L):

        # Call fitting function to get integral and uncertainty
        bkg_fn_fit = ROOT.TF1("bkg_fn_fit", SB_function, d.GetBinLowEdge(1), d.GetBinLowEdge(d.GetNbinsX()) + d.GetBinWidth(1), 3)
        data_fit_result = d.Fit(bkg_fn_fit, "S0")
        cov_data =  data_fit_result.GetCovarianceMatrix()
        

        bkg_fn_fit_Nm1 = ROOT.TF1("bkg_fn_fit_Nm1", SB_function, d.GetBinLowEdge(1), d.GetBinLowEdge(d.GetNbinsX()) + d.GetBinWidth(1), 3)
        data_fit_result_Nm1 = d_Nm1.Fit(bkg_fn_fit_Nm1, "S0")
        cov_data_Nm1 =  data_fit_result_Nm1.GetCovarianceMatrix()

        # For MC closure Call fitting function to get integral and uncertainty
        bkg_fn_fit_mc = ROOT.TF1("bkg_fn_fit_mc", SB_function, d.GetBinLowEdge(1), d.GetBinLowEdge(d.GetNbinsX()) + d.GetBinWidth(1), 3)
        mc.Fit(bkg_fn_fit_mc, "S0")

        # For MC closure Call fitting function to get integral and uncertainty
        bkg_fn_fit_mc_Nm1 = ROOT.TF1("bkg_fn_fit_mc_Nm1", SB_function, mc_Nm1.GetBinLowEdge(1), mc_Nm1.GetBinLowEdge(mc_Nm1.GetNbinsX()) + mc_Nm1.GetBinWidth(1), 3)
        mc_Nm1.Fit(bkg_fn_fit_mc_Nm1, "S0")

        #Background function with full range
        bkg_fn = ROOT.TF1("bkg_fn", "pol2", d.GetBinLowEdge(1), d.GetBinLowEdge(d.GetNbinsX()) + d.GetBinWidth(1))
        bkg_fn_Nm1 = ROOT.TF1("bkg_fn_Nm1", "pol2", d_Nm1.GetBinLowEdge(1), d_Nm1.GetBinLowEdge(d_Nm1.GetNbinsX()) + d_Nm1.GetBinWidth(1))
        bkg_fn_mc = ROOT.TF1("bkg_fn_mc", "pol2", mc.GetBinLowEdge(1), mc.GetBinLowEdge(mc.GetNbinsX()) + mc.GetBinWidth(1))
        bkg_fn_mc_Nm1 = ROOT.TF1("bkg_fn_mc_Nm1", "pol2", mc_Nm1.GetBinLowEdge(1), mc_Nm1.GetBinLowEdge(mc_Nm1.GetNbinsX()) + mc_Nm1.GetBinWidth(1))

        #Pass parameters from ammended range bkg fn to full range bkg fn
        bkg_fn.SetParameters(bkg_fn_fit.GetParameters())
        bkg_fn.SetParErrors(bkg_fn_fit.GetParErrors())
        bkg_fns += [bkg_fn]
        data_bkg_in_peak = ROOT.Double()
        data_bkg_in_peak = bkg_fn.Integral(LB, UB) / d.GetBinWidth(1)
        data_total_in_peak = d.Integral(d.GetXaxis().FindBin(LB) + 1, d.GetXaxis().FindBin(UB) - 1)
        
        #Subtract bkg chi2 dist from peak region by weighting by SB fit to data
        #Fraction of peak region vs whole range
        frac2 = data_bkg_in_peak/(d.Integral(1,d.GetXaxis().FindBin(LB)) + d.Integral(d.GetXaxis().FindBin(UB),-1))
        data_varTest_OutSB = ROOT.TH1F()
        data_varTest_OutSB = data_chi2OutSB_R.Clone()
        data_varTest_OutSB.Add(data_chi2OutSB_L)
        data_chi2InSB.Add(data_varTest_OutSB,-1*frac2)
        print( "Fraction of bkg in peak for data: "+ str(data_bkg_in_peak/data_total_in_peak))

        #Pass parameters from ammended range bkg fn to full range bkg fn for N-1 Denom of efficiency
        bkg_fn_Nm1.SetParameters(bkg_fn_fit_Nm1.GetParameters())
        data_bkg_in_peak_Nm1 = ROOT.Double()
        data_bkg_in_peak_Nm1 = bkg_fn_Nm1.Integral(LB, UB) / d_Nm1.GetBinWidth(1)
        data_total_in_peak_Nm1 = d_Nm1.Integral(d_Nm1.GetXaxis().FindBin(LB) + 1, d_Nm1.GetXaxis().FindBin(UB) - 1)
        
        #Subtract bkg chi2 dist from peak region by weighting by SB fit to data
        #Fraction of peak region vs whole range
        frac2_Nm1 = data_bkg_in_peak_Nm1/(d_Nm1.Integral(1,d_Nm1.GetXaxis().FindBin(LB)) + d_Nm1.Integral(d.GetXaxis().FindBin(UB),-1))
        data_varTest_OutSB = ROOT.TH1F()
        data_varTest_OutSB = data_chi2OutSB_R.Clone()
        data_varTest_OutSB.Add(data_chi2OutSB_L)
        data_chi2InSB.Add(data_varTest_OutSB,-1*frac2)
        print( "Fraction of bkg in peak for data N-1: "+ str(data_bkg_in_peak_Nm1/data_total_in_peak_Nm1))

        #Pass parameters from ammended range bkg fn to full range bkg fn for Matched MC
        bkg_fn_mc.SetParameters(bkg_fn_fit_mc.GetParameters())
        data_bkg_in_peak_mc = ROOT.Double()
        data_bkg_in_peak_mc = bkg_fn_mc.Integral(LB, UB) / mc.GetBinWidth(1)
        data_total_in_peak_mc = mc.Integral(mc.GetXaxis().FindBin(LB) + 1, mc.GetXaxis().FindBin(UB) - 1)

        bkg_fn_mc_Nm1.SetParameters(bkg_fn_fit_mc_Nm1.GetParameters())
        data_bkg_in_peak_mc_Nm1 = ROOT.Double()
        data_bkg_in_peak_mc_Nm1 = bkg_fn_mc_Nm1.Integral(LB, UB) / mc_Nm1.GetBinWidth(1)
        data_total_in_peak_mc_Nm1 = mc_Nm1.Integral(mc_Nm1.GetXaxis().FindBin(LB) + 1, mc_Nm1.GetXaxis().FindBin(UB) - 1)

        mc_err = ROOT.Double()
        mc_count = mc.IntegralAndError(1,mc.GetNbinsX(),mc_err)
        mc_err_Nm1 = ROOT.Double()
        mc_Nm1_count = mc_Nm1.IntegralAndError(1,mc_Nm1.GetNbinsX(),mc_err_Nm1)
        #Efficiency and Uncertainty numbers
        if c is "OS":
            OS_err.append(data_total_in_peak)
            OS_err.append(data_total_in_peak_Nm1)
            OS_mc_err.append(mc_err)
            OS_mc_err.append(mc_err_Nm1)
        elif c is "SS":
            SS_err.append(data_total_in_peak)
            SS_err.append(data_total_in_peak_Nm1)
            SS_mc_err.append(mc_err)
            SS_mc_err.append(mc_err_Nm1)
        elif c is "OS-SS":
            #error data
            num = data_total_in_peak - data_bkg_in_peak
            den = data_total_in_peak_Nm1 - data_bkg_in_peak_Nm1
            data_eff = num/den
            #error on events in data
            num_err = (OS_err[0] + SS_err[0] + (bkg_fn.IntegralError(LB, UB,bkg_fn.GetParameters(),  cov_data.GetMatrixArray()) / d.GetBinWidth(1))**2) ** 0.5
            den_err = (OS_err[1] + SS_err[1] + (bkg_fn_Nm1.IntegralError(LB, UB, bkg_fn_Nm1.GetParameters(),  cov_data_Nm1.GetMatrixArray()) / d_Nm1.GetBinWidth(1))**2) ** 0.5
            data_err = data_eff * (((num_err/num)**2 + (den_err/den)**2)**0.5)
            #mc eff
            #mc_eff = (data_total_in_peak_mc - data_bkg_in_peak_mc)/(data_total_in_peak_mc_Nm1 - data_bkg_in_peak_mc_Nm1)
            mc_eff = mc_count/mc_Nm1_count
            num_mc_err = (OS_mc_err[0]**2 + SS_mc_err[0]**2)**0.5
            den_mc_err = (OS_mc_err[1]**2 + SS_mc_err[1]**2)**0.5
            mc_eff_err = mc_eff * (((num_mc_err/mc_count)**2+(den_mc_err/mc_Nm1_count)**2))**0.5
            #ratio Error
            sf_ratio = data_eff/mc_eff
            sf_err = sf_ratio * ((mc_eff_err/mc_eff)**2+(data_err/data_eff)**2)**0.5
            print(str(data_eff) + " " + str(data_err))
            print("data_eff: " + str(data_eff))
            print("data_err: " + str(data_err))
            print("mc_eff: " + str(mc_eff))
            print("mc_err: " + str(mc_eff_err))
            print("SF: " + str(sf_ratio))
            print("SF_err: " + str(sf_err))
            print("data err: num, den, fit, " + str((OS_err[0] + SS_err[0])**0.5) + " : " + str(bkg_fn.IntegralError(LB, UB, bkg_fn.GetParameters(),  cov_data.GetMatrixArray()) / d.GetBinWidth(1)))
            print(str(d.GetXaxis().FindBin(LB))+" : "+str(d.GetXaxis().FindBin(LB-0.01))+" : "+str(d.GetXaxis().FindBin(LB-0.02)))
             
#        ## Repeat of data but mc for closure
#        #Pass parameters from ammended range bkg fn to full range bkg fn
#        bkg_fn_mc.SetParameters(bkg_fn_fit_mc.GetParameters())
#        mc_bkg_in_peak = ROOT.Double()
#        # Events and Uncertainty in SB region
#        #Integral is the amount of data in Signal region
#        mc_bkg_in_peak = bkg_fn_mc.Integral(LB, UB) / d.GetBinWidth(1)
#        mc_total_in_peak = mc_all.Integral(d.GetXaxis().FindBin(LB) + 1, d.GetXaxis().FindBin(UB) - 1)
#        
#        #Subtract bkg chi2 dist from peak region by weighting by SB fit to data
#        #Fraction of peak region vs whole range
#        frac2_mc = mc_bkg_in_peak/(mc_all.Integral(1,d.GetXaxis().FindBin(LB)) + d.Integral(d.GetXaxis().FindBin(UB),-1))
#        mc_all_varTest_OutSB = ROOT.TH1F()
#        mc_all_varTest_OutSB = mc_all_varTest_outSB_R.Clone()
#        mc_all_varTest_OutSB.Add(mc_all_varTest_outSB_L)
#        mc_all_varTest_inSB.Add(mc_all_varTest_OutSB,-1*frac2_mc)
#        print( "Fraction of bkg in peak for MC: "+ str(mc_bkg_in_peak/mc_total_in_peak))
#        
#
        ### Making plots
        canv_1 = ROOT.TCanvas(c, c, 800, 600)
        #d.GetXaxis().SetRangeUser(135, 180)
        d.GetYaxis().SetRangeUser(d.GetMinimum(), d.GetMaximum() * 1.3)
        d.GetXaxis().SetTitle("D Mass [GeV]")
        l1 = ROOT.TLatex()
        l1.SetTextFont(72)
        l1.SetTextSize(28)
        l1.DrawLatex(0.18, .85, "ATLAS Internal")
        leg1 = ROOT.TLegend(.65,.7,.9,.95)
        leg1.AddEntry(d, "Data")
        leg1.AddEntry(bkg_fn, "SB fit to bkg","l")
        d.Draw("")
        bkg_fn.SetLineColor(2)
        bkg_fn.Draw("same")
        leg1.Draw()
        #l1.DrawLatex(0.18, .85, "ATLAS Internal")
        canv_1.Print(c + "_" + var_test + "_" + "data_bkg_fit" + ntags[0] +".png")

        canv_5 = ROOT.TCanvas(c, c, 800, 600)
        #d_Nm1.GetXaxis().SetRangeUser(135, 180)
        d_Nm1.GetYaxis().SetRangeUser(d_Nm1.GetMinimum(), d_Nm1.GetMaximum() * 1.3)
        d_Nm1.GetXaxis().SetTitle("D Mass [GeV]")
        l1 = ROOT.TLatex()
        l1.SetTextFont(72)
        l1.SetTextSize(28)
        l1.DrawLatex(0.18, .85, "ATLAS Internal")
        leg5 = ROOT.TLegend(.65,.7,.9,.95)
        leg5.AddEntry(d, "Data")
        leg5.AddEntry(bkg_fn, "SB fit to bkg","l")
        d_Nm1.Draw("")
        bkg_fn_Nm1.SetLineColor(2)
        bkg_fn_Nm1.Draw("same")
        leg5.Draw()
        #l1.DrawLatex(0.18, .85, "ATLAS Internal")
        canv_5.Print(c + "_" + var_test + "_" + "_data_bkg_fit_Nm1" + ntags[0] +".png")

        bkgHist = bkg_fn.CreateHistogram()
        bkgHist_Nm1 = bkg_fn_Nm1.CreateHistogram()
        bkgHist_mc = bkg_fn.CreateHistogram()
        bkgHist_mc_Nm1 = bkg_fn_Nm1.CreateHistogram()
        print("bin_low_edge: " + str(bkgHist.GetBinLowEdge(1)) + " : " + str(d.GetBinLowEdge(1)))
        print("Nm1_bin_low_edge: " + str(bkgHist_Nm1.GetBinLowEdge(1)) + " : " + str(d_Nm1.GetBinLowEdge(1)))
        print("bin_N: " + str(bkgHist.GetNbinsX()) + " : " + str(d.GetNbinsX()))
        print("Nm1_bin_N: " + str(bkgHist_Nm1.GetNbinsX()) + " : " + str(d_Nm1.GetNbinsX()))
        
        canv_12 = ROOT.TCanvas(c, c, 800, 600)
        #d.GetXaxis().SetRangeUser(135, 180)
        d.Add(bkgHist,-1)
        d.GetYaxis().SetRangeUser(d.GetMinimum(), d.GetMaximum() * 1.3)
        d.GetXaxis().SetTitle("D Mass [GeV]")
        leg12 = ROOT.TLegend(.65,.7,.9,.95)
        leg12.AddEntry(d, "Data (SB Subtracted)")
        leg12.AddEntry(mc, "Matched MC","l")
        d.Draw("")
        #mc.Add(bkgHist_mc, -1)
        #mc.Scale(d.Integral()/mc.Integral())
        print(str(d.Integral()) + " : " + str(mc.Integral()))
        mc.SetLineColor(2)
        mc.Draw("histsame")
        leg12.Draw()
        #l1.DrawLatex(0.18, .85, "ATLAS Internal")
        canv_12.Print(c + "_" + var_test + "_" + "_data_mc_comp" + ntags[0] +".png")

        canv_52 = ROOT.TCanvas(c, c, 800, 600)
        d_Nm1.Add(bkgHist_Nm1,-1)
        d_Nm1.GetYaxis().SetRangeUser(d_Nm1.GetMinimum(), d_Nm1.GetMaximum() * 1.3)
        d_Nm1.GetXaxis().SetTitle("D Mass [GeV]")
        leg52 = ROOT.TLegend(.65,.7,.9,.95)
        leg52.AddEntry(d, "Data (SB Subtracted)")
        leg52.AddEntry(mc, "Matched MC","l")
        d_Nm1.Draw("")
        mc_Nm1.SetLineColor(2)
        mc_Nm1.Draw("histsame")
        leg52.Draw()
        #l1.DrawLatex(0.18, .85, "ATLAS Internal")
        canv_52.Print(c + "_" + var_test + "_" + "_data_mc_comp_Nm1" + ntags[0] +".png")

#        canv = ROOT.TCanvas(c, c, 800, 600)
#        #if var_test is "Lxy":
#        #    canv.cd(1).SetLogy()
#        chi2_fn = ROOT.TF1("chi2_fn", chi2_function, 0, 20.,0)   
#        data_chi2InSB.Rebin(2)  
#        mc_chi2.Rebin(2)   
#        data_chi2InSB.Scale(1.0/(data_chi2InSB.Integral() * data_chi2InSB.GetBinWidth(1)))
#        mc_chi2.Scale(data_chi2InSB.Integral()/mc_chi2.Integral())
#        mc_chi2.SetLineColor(7)
#        #chi2_fn.SetLineColor(2)
#        data_chi2InSB.GetXaxis().SetTitle(xlabel)
#        data_chi2InSB.GetYaxis().SetTitle("a.u.")
#        #l1 = ROOT.TLatex()
#        #l1.SetTextFont(72)
#        #l1.SetTextSize(28)
#        #l1.DrawLatex(0.18, .85, "ATLAS Internal")
#        leg2 = ROOT.TLegend(.65,.7,.9,.95)
#        leg2.AddEntry(data_chi2InSB, "Data (SB Subtracted)")
#        leg2.AddEntry(mc_chi2, "Matched MC (W+jets & top)","l")
#        #leg2.AddEntry(chi2_fn,"#chi^{2} w/ 3 d.o.f.","l")
#        y_upper_lim = 1.3 * max(data_chi2InSB.GetMaximum(), mc_chi2.GetMaximum())
#        data_chi2InSB.GetYaxis().SetRangeUser(data_chi2InSB.GetMinimum(), y_upper_lim)
#        data_chi2InSB.Draw("")
#        mc_chi2.Draw("histsame")
#        #chi2_fn.Draw("same")
#        leg2.Draw()
#        canv.Print(c + "_data_"+ var_test +".png")

#        
#        #Cumulative chi2 for data
#        data_varTest_cum = ROOT.TH1F()
#        mc_varTest_cum = ROOT.TH1F()
#        data_varTest_cum = data_chi2InSB.GetCumulative(isFwd)
#        data_varTest_cum.Scale(1./data_chi2InSB.Integral())
#        mc_varTest_cum = mc_chi2.GetCumulative(isFwd)
#        mc_varTest_cum.Scale(1./mc_chi2.Integral())
#
#        canv_3 = ROOT.TCanvas(c, c, 800, 600)
#        y_upper_lim = 1.3 * max(data_varTest_cum.GetMaximum(), mc_varTest_cum.GetMaximum())
#        data_varTest_cum.GetYaxis().SetRangeUser(0.0, y_upper_lim)
#        #data_varTest_cum.Divide(mc_varTest_cum)
#        data_varTest_cum.GetYaxis().SetTitle("Cut Efficiency")
#        data_varTest_cum.GetXaxis().SetRangeUser(cut_vals[0],cut_vals[len(cut_vals)-1])
#        data_varTest_cum.Draw("hist")
#        mc_varTest_cum.SetLineColor(2)
#        mc_varTest_cum.Draw("histsame")
#        cut_graph_errors.Draw("psame")
#        leg4 = ROOT.TLegend(.65,.8,.9,.95)
#        leg4.AddEntry(data_varTest_cum, "Data Chi2", "l")
#        leg4.AddEntry(mc_varTest_cum, "Matched MC Chi2","l")
#        leg4.AddEntry(cut_graph_errors, "Cut Scan Efficiency","p")
#        leg4.Draw()
#        canv_3.Print(c + "_data_cum_"+ var_test +".png")
#        #Division hist 
#        data_varTest_cum.Divide(mc_varTest_cum)
#        data_varTest_cum.GetYaxis().SetRangeUser(0.0,data_varTest_cum.GetMaximum())
#        data_varTest_cum.GetYaxis().SetTitle("Cumulative SB Subtraced Data/ Cumulative Matched MC")
#        data_varTest_cum.Draw("hist")
#        canv_3.Print(c + "_data_cum_div_"+ var_test +".png")
#
#
#        ##MC Closure
#        canv_2 = ROOT.TCanvas(c, c, 800, 600)
#        #if var_test is "Lxy":
#        #    canv_2.cd(1).SetLogy()
#        mc_all_varTest_inSB.Rebin(2)
#        #mc_all_varTest_inSB.Scale(chi2_fn.Integral(0.,20.)/(mc_all_varTest_inSB.Integral() * mc_all_varTest_inSB.GetBinWidth(1)))
#        mc_chi2.Scale(mc_all_varTest_inSB.Integral()/mc_chi2.Integral())
#        data_chi2InSB.GetXaxis().SetTitle(xlabel)
#        data_chi2InSB.GetYaxis().SetTitle("a.u.")
#        #l1 = ROOT.TLatex()
#        #l1.SetTextFont(72)
#        #l1.SetTextSize(28)
#        #l1.DrawLatex(0.18, .85, "ATLAS Internal")
#        leg3 = ROOT.TLegend(.65,.7,.9,.95)
#        leg3.AddEntry( mc_all_varTest_inSB, "MC (SB Subtracted)")
#        leg3.AddEntry(mc_chi2, "Matched MC (W+jets & top)","l")
#        #leg3.AddEntry(chi2_fn,"#chi^{2} w/ 3 d.o.f.","l")
#        y_upper_lim = 1.3 * max(mc_all_varTest_inSB.GetMaximum(), mc_chi2.GetMaximum())
#        mc_all_varTest_inSB.GetYaxis().SetRangeUser(mc_all_varTest_inSB.GetMinimum(), y_upper_lim)
#        mc_all_varTest_inSB.GetXaxis().SetTitle(xlabel)
#        mc_all_varTest_inSB.Draw("")
#        mc_chi2.Draw("histsame")
#        #chi2_fn.SetLineColor(2)
#        #chi2_fn.Draw("same")
#        leg3.Draw()
#        canv_2.Print(c + "_mc_"+ var_test +".png")
#        print(data_chi2InSB.Integral())
#        #print(chi2_fn.Integral(0.,20.))
#