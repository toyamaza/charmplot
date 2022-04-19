#!/usr/bin/env python
import os
import ROOT

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasLabels.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasUtils.C"))
ROOT.SetAtlasStyle()

# Include custom PDFs
ROOT.gROOT.LoadMacro(os.path.join(os.path.dirname(__file__), "RooTwoSidedCBShape.cc"))

# Make output folder
if not os.path.isdir("fits"):
    os.makedirs("fits")

# meson_dict = {
#   'DstarPlus' = '61'
#   'DstarMinus' = '71'
#   'DPlus' = '64'
#   'DMinus' = '74'
# }

# spg dict map -> systematic : [root_file, distribution_name, inclusive sigma, pt_bin1 sigma, ...]
spg_dict = {
    "DstarPlus_nominal": ['test999961_01_00_01_FTFP', 'Dstar_pTVsDeltaM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DstarPlus_TRACK_EFF_Overal": ['test999961_01_00_02_FTFP', 'Dstar_pTVsDeltaM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DstarPlus_TRACK_EFF_IBL": ['test999961_01_00_03_FTFP', 'Dstar_pTVsDeltaM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DstarPlus_TRACK_EFF_PP0": ['test999961_01_00_04_FTFP', 'Dstar_pTVsDeltaM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DstarPlus_TRACK_EFF_QGSP": ['test999961_01_00_01_QGSP', 'Dstar_pTVsDeltaM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DstarMinus_nominal": ['test999971_01_00_01_FTFP', 'Dstar_pTVsDeltaM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DstarMinus_TRACK_EFF_Overal": ['test999971_01_00_02_FTFP', 'Dstar_pTVsDeltaM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DstarMinus_TRACK_EFF_IBL": ['test999971_01_00_03_FTFP', 'Dstar_pTVsDeltaM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DstarMinus_TRACK_EFF_PP0": ['test999971_01_00_04_FTFP', 'Dstar_pTVsDeltaM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DstarMinus_TRACK_EFF_QGSP": ['test999971_01_00_01_QGSP', 'Dstar_pTVsDeltaM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DPlus_nominal": ['test999964_01_00_01_FTFP','Dplus_pTVsM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DPlus_TRACK_EFF_Overal": ['test999964_01_00_02_FTFP', 'Dplus_pTVsM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DPlus_TRACK_EFF_IBL": ['test999964_01_00_03_FTFP', 'Dplus_pTVsM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DPlus_TRACK_EFF_PP0": ['test999964_01_00_04_FTFP', 'Dplus_pTVsM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DPlus_TRACK_EFF_QGSP": ['test999964_01_00_01_QGSP', 'Dplus_pTVsM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DMinus_nominal": ['test999974_01_00_01_FTFP', 'Dplus_pTVsM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DMinus_TRACK_EFF_Overal": ['test999974_01_00_02_FTFP', 'Dplus_pTVsM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DMinus_TRACK_EFF_IBL": ['test999974_01_00_03_FTFP', 'Dplus_pTVsM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DMinus_TRACK_EFF_PP0": ['test999974_01_00_04_FTFP', 'Dplus_pTVsM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    "DMinus_TRACK_EFF_QGSP": ['test999974_01_00_01_QGSP', 'Dplus_pTVsM', 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
}

ptbins = ["inclusive", "pt_bin1", "pt_bin2", "pt_bin3", "pt_bin4", "pt_bin5"]


def main():

    # output file
    out = ROOT.TFile("fits/spg_fit_track_eff_sys.root", "RECREATE")

    dict_pt_index = 0

    # loop
    for sys in spg_dict.keys():

        for ptbin in ptbins:

            print(f'\n\nXXX sys key = {sys} XXX\n\n')
            print(f'\n\nXXX ptbin = {ptbin} XXX\n\n')

            f = ROOT.TFile(f'{spg_dict[sys][0]}.root', "READ")

            # Get the histogram
            f.cd()
            h_tmp = f.Get(spg_dict[sys][1])
            print(h_tmp.GetName())
            h_2D_tmp = h_tmp.Clone(f"{h_tmp.GetName()}_2D_clone")

            # pT bins [8, 12, 20, 40, 80, 150]
            if ptbin == "inclusive":
                h = h_2D_tmp.ProjectionX(f"{h_tmp.GetName()}_{ptbin}_clone")
                dict_pt_index = 2
            elif ptbin == "pt_bin1":
                h = h_2D_tmp.ProjectionX(f"{h_tmp.GetName()}_{ptbin}_clone", 9, 12)
                dict_pt_index = 3
            elif ptbin == "pt_bin2":
                h = h_2D_tmp.ProjectionX(f"{h_tmp.GetName()}_{ptbin}_clone", 13, 20)
                dict_pt_index = 4
            elif ptbin == "pt_bin3":
                h = h_2D_tmp.ProjectionX(f"{h_tmp.GetName()}_{ptbin}_clone", 21, 40)
                dict_pt_index = 5
            elif ptbin == "pt_bin4":
                h = h_2D_tmp.ProjectionX(f"{h_tmp.GetName()}_{ptbin}_clone", 41, 80)
                dict_pt_index = 6
            elif ptbin == "pt_bin5":
                h = h_2D_tmp.ProjectionX(f"{h_tmp.GetName()}_{ptbin}_clone", 81, 150)
                dict_pt_index = 7

            if "Dstar" in sys:
                # Observable, i.e. the invariant mass
                x = ROOT.RooRealVar("x", "m(D*^{+}-D^{0})", 140, 150, "MeV")
                x_arg_list = ROOT.RooArgList(x)

                # Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
                dh = ROOT.RooDataHist("dh", "dh", x_arg_list, ROOT.RooFit.Import(h))

                # Signal
                if ptbin != "pt_bin5" and ptbin != "pt_bin4":
                    aLo = ROOT.RooRealVar("aLo", "aLo", 1.375, 1.375, 1.375)
                    aHi = ROOT.RooRealVar("aHi", "aHi", 1.4, 1.4, 1.4)
                elif ptbin == "pt_bin4":
                    aLo = ROOT.RooRealVar("aLo", "aLo", 1.375, 1.375, 1.375)
                    aHi = ROOT.RooRealVar("aHi", "aHi", 1.375, 1.375, 1.375)
                elif ptbin == "pt_bin5":
                    aLo = ROOT.RooRealVar("aLo", "aLo", 1.375, 1.375, 1.375)
                    aHi = ROOT.RooRealVar("aHi", "aHi", 1.1, 1.1, 1.1)
                mean = ROOT.RooRealVar("mean", "mean", 145.5, 140.0, 150.0)
                # aLo = ROOT.RooRealVar("aLo", "aLo", 1, 0, 10)
                # aHi = ROOT.RooRealVar("aHi", "aHi", 1, 0, 10)
                sigma = ROOT.RooRealVar("sigma", "sigma", 0.65, 0., 2.0)
                nLo = ROOT.RooRealVar("nLo", "nLo", 99, 0, 100)
                nHi = ROOT.RooRealVar("nHi", "nHi", 99, 0, 100)
                sig = ROOT.RooTwoSidedCBShape("signal", "signal", x, mean, sigma, aLo, nLo, aHi, nHi)
            else:
                # Observable, i.e. the invariant mass
                x = ROOT.RooRealVar("x", "m(D^{+})", 1.7, 2.05, "GeV")
                x_arg_list = ROOT.RooArgList(x)

                # Create a binned dataset that imports contents of TH1 and associates its contents to observable 'x'
                dh = ROOT.RooDataHist("dh", "dh", x_arg_list, ROOT.RooFit.Import(h))

                # Signal
                if ptbin != "pt_bin5" and ptbin != "pt_bin1":
                    mean = ROOT.RooRealVar("mean", "mean", 1.870, 1.86, 1.88)
                    aLo = ROOT.RooRealVar("aLo", "aLo", 1.25, 1.25, 1.25)
                    aHi = ROOT.RooRealVar("aHi", "aHi", 1.25, 1.25, 1.25)
                elif ptbin == "pt_bin1":
                    mean = ROOT.RooRealVar("mean", "mean", 1.870, 1.86, 1.88)
                    aLo = ROOT.RooRealVar("aLo", "aLo", 1.21, 1.21, 1.21)
                    aHi = ROOT.RooRealVar("aHi", "aHi", 1.21, 1.21, 1.21)
                elif ptbin == "pt_bin5":
                    mean = ROOT.RooRealVar("mean", "mean", 1.873, 1.86, 1.88)
                    aLo = ROOT.RooRealVar("aLo", "aLo", 1.4, 1.4, 1.4)
                    aHi = ROOT.RooRealVar("aHi", "aHi", 1.4, 1.4, 1.4)
                sigma = ROOT.RooRealVar("sigma", "sigma", 0.021, 0.01, 0.03)
                # nLo = ROOT.RooRealVar("nLo", "nLo", 4, 1, 10)
                # nHi = ROOT.RooRealVar("nHi", "nHi", 4, 1, 10)
                nLo = ROOT.RooRealVar("nLo", "nLo", 99, 0, 100)
                nHi = ROOT.RooRealVar("nHi", "nHi", 99, 0, 100)
                sig = ROOT.RooTwoSidedCBShape("signal", "signal", x, mean, sigma, aLo, nLo, aHi, nHi)

            # Do the fit
            _ = sig.fitTo(dh, ROOT.RooFit.Minimizer("Minuit2", "Migrad"),
                          ROOT.RooFit.SumW2Error(ROOT.kTRUE),
                          ROOT.RooFit.Save(),
                          ROOT.RooFit.PrintLevel(1))

            s_component = ROOT.RooArgSet(sig)
            if "Dstar" in sys:
                frame = x.frame(ROOT.RooFit.Title("Fit D* mass"))
            else:
                frame = x.frame(ROOT.RooFit.Title("Fit D+ mass"))

            # data
            dh.plotOn(frame, ROOT.RooFit.MarkerColor(ROOT.kBlack), ROOT.RooFit.Name("data"), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))

            # signal
            sig.plotOn(frame, ROOT.RooFit.Components(s_component), ROOT.RooFit.LineColor(
                ROOT.kRed), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.Name("s"))

            # canvas
            c = ROOT.TCanvas(sys, sys, 1000, 1000)
            frame.Draw()

            # print text
            ROOT.ATLASLabel(0.17, 0.90, "Internal", 1)
            ROOT.myText(0.17, 0.9 - 1 * 0.06, 1, "#sqrt{s} = 13 TeV")
            ROOT.myText(0.17, 0.9 - 2 * 0.06, 1, f"{sys.split('_')[0]}")
            ROOT.myText(0.17, 0.9 - 3 * 0.06, 1, f"{sys.split('_')[-1]}")
            ROOT.myText(0.17, 0.9 - 4 * 0.06, 1, f"{ptbin}")
            ROOT.myText(0.68, 0.9 - 0 * 0.06, 1, f"aHi: {aHi.getVal():.4f}")
            ROOT.myText(0.68, 0.9 - 1 * 0.06, 1, f"aLo: {aLo.getVal():.4f}")
            ROOT.myText(0.68, 0.9 - 2 * 0.06, 1, f"mean: {mean.getVal():.4f}")
            ROOT.myText(0.68, 0.9 - 3 * 0.06, 1, f"nHi: {nHi.getVal():.4f}")
            ROOT.myText(0.68, 0.9 - 4 * 0.06, 1, f"nLo: {nLo.getVal():.4f}")
            if "Dstar" in sys:
                ROOT.myText(0.68, 0.9 - 5 * 0.06, 1, f"sigma: {sigma.getVal():.6f}")
            else:
                ROOT.myText(0.68, 0.9 - 5 * 0.06, 1, f"sigma: {1000*sigma.getVal():.6f}")

            # save to output
            h_out = ROOT.TH1D(f"pars_{sys}_{ptbin}", f"pars_{sys}_{ptbin}", 6, -0.5, 5.5)
            h_out.SetBinContent(1, aHi.getVal())
            h_out.SetBinContent(2, aLo.getVal())
            h_out.SetBinContent(3, mean.getVal())
            h_out.SetBinContent(4, nHi.getVal())
            h_out.SetBinContent(5, nLo.getVal())
            h_out.SetBinContent(6, sigma.getVal())
            h_out.SetBinError(1, aHi.getError())
            h_out.SetBinError(2, aLo.getError())
            h_out.SetBinError(3, mean.getError())
            h_out.SetBinError(4, nHi.getError())
            h_out.SetBinError(5, nLo.getError())
            h_out.SetBinError(6, sigma.getError())
            h_out.GetXaxis().SetBinLabel(1, aHi.GetName())
            h_out.GetXaxis().SetBinLabel(2, aLo.GetName())
            h_out.GetXaxis().SetBinLabel(3, mean.GetName())
            h_out.GetXaxis().SetBinLabel(4, nHi.GetName())
            h_out.GetXaxis().SetBinLabel(5, nLo.GetName())
            h_out.GetXaxis().SetBinLabel(6, sigma.GetName())
            out.cd()
            h_out.Write()

            # Save the fit into pdf/png files
            c.Print(f"fits/{sys}_{ptbin}.pdf")
            c.Print(f"fits/{sys}_{ptbin}.png")

            # Save sigma value
            spg_dict[sys][dict_pt_index] = sigma.getVal()

            f.Close()
        # ptbin end
    # sys end

    # inclusive_sys - sys : [inclusive sigma, pt_bin1 sigma, ...]
    inclusive_sys = {
        "Dstar_nominal": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "Dstar_TRACK_EFF_Overal": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "Dstar_TRACK_EFF_IBL": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "Dstar_TRACK_EFF_PP0": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "Dstar_TRACK_EFF_QGSP": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "Dplus_nominal": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "Dplus_TRACK_EFF_Overal": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "Dplus_TRACK_EFF_IBL": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "Dplus_TRACK_EFF_PP0": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        "Dplus_TRACK_EFF_QGSP": [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    }

    sigma_diff = 0

    # Calculate sigma diff for systematic variations
    for sys in spg_dict.keys():
        print("\n")
        print(sys)
        for i in range(len(spg_dict[sys])):
            # Skip nominal fits and first two entries of dictionary
            if i < 2:
                continue
            if "nominal" in sys:
                if "Dstar" in sys:
                    inclusive_sys["Dstar_nominal"][i - 2] += 0.5 * spg_dict[sys][i]
                elif "DPlus" in sys:
                    inclusive_sys["Dplus_nominal"][i - 2] += 0.5 * 1000 * spg_dict[sys][i]
                elif "DMinus" in sys:
                    inclusive_sys["Dplus_nominal"][i - 2] += 0.5 * 1000 * spg_dict[sys][i]
                continue
            print(ptbins[i - 2])

            # Get name of systematic for inclusive sys dictionary
            if "Overal" in sys:
                i_sys = "TRACK_EFF_Overal"
            elif "IBL" in sys:
                i_sys = "TRACK_EFF_IBL"
            elif "PP0" in sys:
                i_sys = "TRACK_EFF_PP0"
            elif "QGSP" in sys:
                i_sys = "TRACK_EFF_QGSP"

            if "DstarPlus" in sys:
                sigma_diff = ROOT.TMath.Sqrt(abs(spg_dict[sys][i]**2 - spg_dict["DstarPlus_nominal"][i]**2))
                # sigma_diff_var = ROOT.RooRealVar(f"{sys}_{ptbins[i - 2]}_sigma_diff", f"{sys}_{ptbins[i - 2]}_sigma_diff", sigma_diff)
                # sigma_diff_var.Write()
                inclusive_sys["Dstar_" + i_sys][i - 2] += 0.5 * spg_dict[sys][i]
                print(f'nominal sigma = {spg_dict["DstarPlus_nominal"][i]}')
                print(f'sys sigma = {spg_dict[sys][i]}')
                print(f'sigma diff = {sigma_diff}')
            elif "DstarMinus" in sys:
                sigma_diff = ROOT.TMath.Sqrt(abs(spg_dict[sys][i]**2 - spg_dict["DstarMinus_nominal"][i]**2))
                # sigma_diff_var = ROOT.RooRealVar(f"{sys}_{ptbins[i - 2]}_sigma_diff", f"{sys}_{ptbins[i - 2]}_sigma_diff", sigma_diff)
                # sigma_diff_var.Write()
                inclusive_sys["Dstar_" + i_sys][i - 2] += 0.5 * spg_dict[sys][i]
                print(f'nominal sigma = {spg_dict["DstarMinus_nominal"][i]}')
                print(f'sys sigma = {spg_dict[sys][i]}')
                print(f'sigma diff = {sigma_diff}')
            elif "DPlus" in sys:
                sigma_diff = ROOT.TMath.Sqrt(abs((spg_dict[sys][i] * 1000)**2 - (spg_dict["DPlus_nominal"][i] * 1000)**2))
                # sigma_diff_var = ROOT.RooRealVar(f"{sys}_{ptbins[i - 2]}_sigma_diff", f"{sys}_{ptbins[i - 2]}_sigma_diff", sigma_diff)
                # sigma_diff_var.Write()
                inclusive_sys["Dplus_" + i_sys][i - 2] += 0.5 * 1000 * spg_dict[sys][i]
                print(f'nominal sigma = {spg_dict["DPlus_nominal"][i] * 1000}')
                print(f'sys sigma = {spg_dict[sys][i] * 1000}')
                print(f'sigma diff = {sigma_diff}')
            elif "DMinus" in sys:
                sigma_diff = ROOT.TMath.Sqrt(abs((spg_dict[sys][i] * 1000)**2 - (spg_dict["DMinus_nominal"][i] * 1000)**2))
                # sigma_diff_var = ROOT.RooRealVar(f"{sys}_{ptbins[i - 2]}_sigma_diff", f"{sys}_{ptbins[i - 2]}_sigma_diff", sigma_diff)
                # sigma_diff_var.Write()
                inclusive_sys["Dplus_" + i_sys][i - 2] += 0.5 * 1000 * spg_dict[sys][i]
                print(f'nominal sigma = {spg_dict["DMinus_nominal"][i] * 1000}')
                print(f'sys sigma = {spg_dict[sys][i] * 1000}')
                print(f'sigma diff = {sigma_diff}')
    
    i_sigma_diff = 0
    dstar_tot_sigma_diff = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    dplus_tot_sigma_diff = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    for sys in inclusive_sys.keys():
        print("\n")
        print(sys)
        for i in range(len(inclusive_sys[sys])):
            print(ptbins[i])
            if "nominal" in sys:
                print(f"sigma = {inclusive_sys[sys][i]}")
                continue
            if "Dstar" in sys:
                print(f"sigma = {inclusive_sys[sys][i]}")
                i_sigma_diff = ROOT.TMath.Sqrt(abs((inclusive_sys[sys][i])**2 - (inclusive_sys["Dstar_nominal"][i])**2))
                dstar_tot_sigma_diff[i] += i_sigma_diff**2
                print(f"sigma diff = {i_sigma_diff}")
            elif "Dplus" in sys:
                print(f"sigma = {inclusive_sys[sys][i]}")
                i_sigma_diff = ROOT.TMath.Sqrt(abs((inclusive_sys[sys][i])**2 - (inclusive_sys["Dplus_nominal"][i])**2))
                dplus_tot_sigma_diff[i] += i_sigma_diff**2
                print(f"sigma diff = {i_sigma_diff}")

    dstar_tot_sigma_diff = [ROOT.TMath.Sqrt(diff) for diff in dstar_tot_sigma_diff]
    dplus_tot_sigma_diff = [ROOT.TMath.Sqrt(diff) for diff in dplus_tot_sigma_diff]

    print("\n")
    print(f"dstar total sigma diff = {dstar_tot_sigma_diff}")
    print(f"dplus total sigma diff = {dplus_tot_sigma_diff}")

    for i in range(len(dstar_tot_sigma_diff)):
        sigma_diff_var = ROOT.RooRealVar(f"Dstar_TRACK_EFF_TOT_{ptbins[i]}_sigma_diff", f"Dstar_TRACK_EFF_TOT_{ptbins[i]}_sigma_diff", dstar_tot_sigma_diff[i])
        sigma_diff_var.Write()

    for i in range(len(dplus_tot_sigma_diff)):
        sigma_diff_var = ROOT.RooRealVar(f"Dplus_TRACK_EFF_TOT_{ptbins[i]}_sigma_diff", f"Dplus_TRACK_EFF_TOT_{ptbins[i]}_sigma_diff", dplus_tot_sigma_diff[i])
        sigma_diff_var.Write()

    out.Close()


if __name__ == "__main__":

    # run
    main()
