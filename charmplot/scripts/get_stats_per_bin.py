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

if __name__ == "__main__":

    file = ROOT.TFile("histograms.root")

    mc_samples = ["Multijet_MatrixMethod","Wjets_emu_Rest", "Wjets_emu_Matched", "Top_Rest", "Top_Matched", "Zjets_emu", "Other"]
#    mc_samples = ["Wjets_emu_Charm", "Wjets_emu_Matched", "Top_Charm", "Top_Matched", "Zjets_emu", "Other"]
#    mc_samples =   ["Wjets_emu_HardMisMatched", "Wjets_emu_Matched","Wjets_emu_NoMatch","Wjets_emu_Charm", "Top_Matched", "Top_Charm", "Top_HardMisMatched" , "Zjets_emu", "Other"]

    data_samples = ["Data"]

    charges = ["OS", "SS", "OS-SS"]
#    charges = ["OS-SS"]
    leps = ["mu", "el"]

    for i in range(0,1):
        var_in = ""
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
                    if first_hist_mc:
                        mc_sum.Clear()
                        mc_sum = file.Get(s + "_" + c + "_2018_" + l + "_SR_Dplus_Dmeson_pt" + var_in).Clone()
                        first_hist_mc = False
                    else:
                        hist = ROOT.TH1F()
                        #print(s + "_" + c + "_2018_" + l + "_SR_Dplus_Dmeson_pt" + var_in)
                        hist = file.Get(s + "_" + c + "_2018_" + l + "_SR_Dplus_Dmeson_pt" + var_in).Clone()
                        #print(s + "_" + c + "_2018_" + l + "_SR_Dplus_Dmeson_pt" + var_in)
                        mc_sum.Add(hist)
                #print()
                for s in data_samples:

                    if first_hist_data:
                        data_sum.Clear()
                        data_sum = file.Get(s + "_" + c + "_2018_" + l + "_SR_Dplus_Dmeson_pt" + var_in).Clone()
                        #print(s + "_" + c + "_2018_" + l + "_SR_Dplus_Dmeson_pt")
                        first_hist_data = False
                    else:
                        hist = file.Get(s + "_" + c + "_2018_" + l + "_SR_Dplus_Dmeson_pt" + var_in).Clone()
                        #print(s + "_" + c + "_2018_" + l + "_SR_Dplus_Dmeson_pt")
                        #print (hist)
                        data_sum.Add(hist)
                #print()
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
                #print("lep: " + l)
                if first_hist_sig:
                    sig_sum.Clear()
                    sig_sum = file.Get("Wjets_emu_Matched_" + c + "_2018_" + l + "_SR_Dplus_Dmeson_pt" + var_in).Clone()
                    first_hist_sig = False
                    #print("el MC: " + str(sig_sum.Integral()))
                else:
                    hist = file.Get("Wjets_emu_Matched_" + c + "_2018_" + l + "_SR_Dplus_Dmeson_pt" + var_in).Clone()
                    sig_sum.Add(hist)
                    #print("mu MC: " + str(hist.Integral()))
                if first_hist_sig_data:
                    sig_sum_data.Clear()
                    sig_sum_data = file.Get("Data_" + c + "_2018_" + l + "_SR_Dplus_Dmeson_pt" + var_in).Clone()
                    first_hist_sig_data = False
                    #print("el MC: " + str(sig_sum_data.Integral()))
                else:
                    hist_data = ROOT.TH1F()
                    hist_data = file.Get("Data_" + c + "_2018_" + l + "_SR_Dplus_Dmeson_pt" + var_in).Clone()
                    sig_sum_data.Add(hist_data)
                    #print("mu MC: " + str(hist_data.Integral()))
            first_hist_sig = True
            first_hist_sig_data = True
            sig_sums += [sig_sum]
            sig_sums_data += [sig_sum_data]

        #print("lep sum: " + str(sig_sums[2].Integral()))
        #print("lep sum: " + str(sig_sums_data[2].Integral()))

        for c, d, mc, sig, sig_data in zip(charges, data_sums, mc_sums, sig_sums, sig_sums_data):

            #var1_bins = [5, 7, 11, 16, 22, 30 ,45, 151]
            var1_bins = [5, 10, 15, 20, 25, 35 ,50, 151]
            root_array1 = array('d',var1_bins)

            pt_var1_hist = ROOT.TH1F("var1","var1",len(var1_bins)-1, root_array1)
            pt_var1_mc_hist = ROOT.TH1F("var1","var1",len(var1_bins)-1, root_array1)
            for i in range(len(var1_bins) - 1):
                pt_var1_hist.SetBinContent(i+1, d.Integral(var1_bins[i],var1_bins[i+1])* 150/58.5)
                pt_var1_mc_hist.SetBinContent(i+1, mc.Integral(var1_bins[i],var1_bins[i+1])* 150/58.5)
                print(str(round(d.Integral(var1_bins[i]+1,var1_bins[i+1]) * 150/58.5)))
                #print("pt [GeV]: " + str(var1_bins[i]) + " to " + str(var1_bins[i+1]) + " : "+ str(round(d.Integral(var1_bins[i],var1_bins[i+1]) * 150/58.5)))

            canv = ROOT.TCanvas(c, c, 800, 600)
            l1 = ROOT.TLatex()
            l1.SetTextFont(72)
            l1.SetTextSize(28)
            l1.DrawLatex(0.18, .85, "ATLAS Internal")
            pt_var1_hist.GetYaxis().SetRangeUser(.5 * min(pt_var1_hist.GetMinimum(),pt_var1_mc_hist.GetMinimum()),1.2*max(pt_var1_hist.GetMaximum(),pt_var1_mc_hist.GetMaximum()))
            pt_var1_hist.Draw("")
            pt_var1_mc_hist.SetLineColor(2)
            pt_var1_mc_hist.Draw("histsame")            
            canv.Print("var1_" + c+".pdf")

