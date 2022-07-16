#!/usr/bin/env python
import ROOT
import os

PREDICTIONS_DIR = "/global/cfs/cdirs/atlas/shapiro/13tev_wd_rivet/forMiha"
PRIORS_DIR = "/global/cfs/cdirs/atlas/wcharm/charmplot_output/Dmeson_2022_06_15/"

THEORY = [
    "ABMP16_5_nnlo",
    "ATLASpdf21_T3",
    "CT18ANNLO",
    "CT18NNLO",
    "MSHT20nnlo_as118",
    "PDF4LHC21_40",
    "NNPDF30_nnlo_as_0118",
    "NNPDF31_nnlo_as_0118",
    "NNPDF40_nnlo_as_01180",
]


def main(options):

    # save chi2
    chi2_dict = {}
    nll_dict = {}

    # quickFit command list
    outfile = open("task_list_NLL_test.txt", "w")

    # loop over stuff
    for var in options.vars.split(","):

        # normalization factors basis
        if options.basis == "default":
            NFS = ["mu_Rc", "mu_Wminus_rel_1", "mu_Wminus_rel_2", "mu_Wminus_rel_3", "mu_Wminus_rel_4",
                   "mu_Wminus_tot", "mu_Wplus_rel_1", "mu_Wplus_rel_2", "mu_Wplus_rel_3", "mu_Wplus_rel_4"]
            suffix = ""
        elif options.basis == "absolute":
            NFS = ["mu_Wminus_1", "mu_Wminus_2", "mu_Wminus_3", "mu_Wminus_4", "mu_Wminus_5",
                   "mu_Wplus_1", "mu_Wplus_2", "mu_Wplus_3", "mu_Wplus_4", "mu_Wplus_5"]
            suffix = "alt_ "

        # cross section priors
        f_prior = ROOT.TFile(os.path.join(PRIORS_DIR, f"fid_eff_{var}_{options.decay.lower()}", "unfolding.root"), "READ")
        h_minus = f_prior.Get(f"Sherpa2211_WplusD_OS-SS_lep_minus_{options.decay}_Kpipi_truth_differential_{var}")
        h_plus = f_prior.Get(f"Sherpa2211_WplusD_OS-SS_lep_plus_{options.decay}_Kpipi_truth_differential_{var}")

        # Multiply scale factor by 0.677 for D*
        br = 1.0
        if options.decay == "Dstar":
            br = 0.677

        sf = 1.0 * br
        if var == "eta":
            sf = 0.5 * br
        priors = {
            "Wminus": h_minus.Integral() / (2. * sf),
            "Wminus_1": h_minus.GetBinContent(1) / (2. * sf),
            "Wminus_2": h_minus.GetBinContent(2) / (2. * sf),
            "Wminus_3": h_minus.GetBinContent(3) / (2. * sf),
            "Wminus_4": h_minus.GetBinContent(4) / (2. * sf),
            "Wminus_5": h_minus.GetBinContent(5) / (2. * sf),
            "Wplus": h_plus.Integral() / (2. * sf),
            "Wplus_1": h_plus.GetBinContent(1) / (2. * sf),
            "Wplus_2": h_plus.GetBinContent(2) / (2. * sf),
            "Wplus_3": h_plus.GetBinContent(3) / (2. * sf),
            "Wplus_4": h_plus.GetBinContent(4) / (2. * sf),
            "Wplus_5": h_plus.GetBinContent(5) / (2. * sf),
        }
        print(f"============ cross section priors for {var} ============")
        for key, val in priors.items():
            print(f"{key}: {val}")

        # normalization factors
        if options.decay == "Dplus":
            f_result = ROOT.TFile(os.path.join(options.input, f"WCharm_lep_obs_OSSS_complete_{suffix}{var}", "Fits", f"WCharm_lep_obs_OSSS_complete_{suffix}{var}.root"), "READ")
        else:
            f_result = ROOT.TFile(os.path.join(options.input, f"WCharm_lep_obs_OSSS_complete_{suffix}{var}", "Fits", f"WCharm_lep_obs_OSSS_complete{suffix}.root"), "READ")
        fr = f_result.Get("nll_simPdf_newasimovData_with_constr")
        nfs = {nf: fr.floatParsFinal().find(nf).getVal() for nf in NFS}
        print("\n============ normalization factors ============")
        for x in NFS:
            print(f"{x:15s} {nfs[x]:.4f} +{fr.floatParsFinal().find(x).getAsymErrorHi():.4f} {fr.floatParsFinal().find(x).getAsymErrorLo():.4f}")

        # covariance matrix
        cov = ROOT.TMatrixD(len(NFS), len(NFS))
        for i in range(len(NFS)):
            for j in range(len(NFS)):
                corr = fr.correlation(NFS[i], NFS[j])
                sigma_i = fr.floatParsFinal().find(NFS[i]).getError()
                sigma_j = fr.floatParsFinal().find(NFS[j]).getError()
                cov[i][j] = sigma_i * sigma_j * corr

        print(f"\n============ Experimental cov for {var} ============")
        cov.Print()

        # quickFit template
        if options.decay == "Dplus":
            quickFit = f"quickFit --minStrat 1 --minTolerance 1e-5 --savefitresult 0 --hesse 0 --saveNP 0 -f {options.input}/WCharm_lep_obs_OSSS_complete_{suffix}{var}/RooStats/WCharm_lep_obs_OSSS_complete_{suffix}{var}_combined_WCharm_lep_obs_OSSS_complete_{suffix}{var}_model.root -w combined -d obsData -p "
        else:
            quickFit = f"quickFit --minStrat 1 --minTolerance 1e-5 --savefitresult 0 --hesse 0 --saveNP 0 -f {options.input}/WCharm_lep_obs_OSSS_complete_{suffix}{var}/RooStats/WCharm_lep_obs_OSSS_complete{suffix}_combined_WCharm_lep_obs_OSSS_complete{suffix}_model.root -w combined -d obsData -p "
        jobName = f"WCharm_lep_obs_OSSS_complete_{suffix}{var}_unconditional"
        for x in NFS:
            quickFit += f"{x},"
        quickFit += f"mu_Top -o {jobName}.root"
        outfile.write(f"{{\"job_name\": \"{jobName}\", \"cmd\": \"{quickFit}\"}}\n")

        # default NLL
        if options.quickfit_input:
            f = ROOT.TFile(os.path.join(options.input, options.quickfit_input, jobName, f"{jobName}.root"), "READ")
            t = f.Get("nllscan")
            for e in t:
                nll_default = e.nll
                break

        # get prediction
        chi2_dict[var] = {}
        nll_dict[var] = {}
        for prediction in THEORY:
            name = "lep_eta" if var == "eta" else "D_pt"
            lep = "minus"
            meson_charge = "plus"
            f_theory = {}
            for lep, meson_charge in zip(["minus", "plus"], ["plus", "minus"]):
                if options.decay == "Dplus":
                    f_name = f"TheoryPredictions_{name}_W{lep}D{meson_charge}.root"
                else:
                    f_name = f"TheoryPredictions_{name}_W{lep}Dstar{meson_charge}.root"
                f = ROOT.TFile(os.path.join(PREDICTIONS_DIR, f_name))
                f_theory[lep] = f
                assert f

            # quickFit template
            if options.decay == "Dplus":
                quickFit = f"quickFit --minStrat 1 --minTolerance 1e-5 --savefitresult 0 --hesse 0 --saveNP 0 -f {options.input}/WCharm_lep_obs_OSSS_complete_{suffix}{var}/RooStats/WCharm_lep_obs_OSSS_complete_{suffix}{var}_combined_WCharm_lep_obs_OSSS_complete_{suffix}{var}_model.root -w combined -d obsData "
            else:
                quickFit = f"quickFit --minStrat 1 --minTolerance 1e-5 --savefitresult 0 --hesse 0 --saveNP 0 -f {options.input}/WCharm_lep_obs_OSSS_complete_{suffix}{var}/RooStats/WCharm_lep_obs_OSSS_complete{suffix}_combined_WCharm_lep_obs_OSSS_complete{suffix}_model.root -w combined -d obsData "

            # read theory graphs
            gr_theory = {}
            gr_theory_rel = {}
            for lep, meson_charge in zip(["minus", "plus"], ["plus", "minus"]):
                gr = f_theory[lep].Get(f"mu_{lep}_{options.decay}_{name}_{prediction}_cross")
                gr_rel = f_theory[lep].Get(f"mu_{lep}_{options.decay}_{name}_{prediction}_norm")
                assert gr, f"mu_{lep}_{options.decay}_{name}_{prediction}_cross"
                gr_theory[lep] = gr
                gr_theory_rel[lep] = gr_rel
            sum_minus = sum([gr_theory['minus'].GetY()[i] for i in range(gr_theory['minus'].GetN())])
            sum_plus = sum([gr_theory['plus'].GetY()[i] for i in range(gr_theory['plus'].GetN())])

            # normalization factor values for the prediction
            prediction_nfs = {}
            if options.basis == "default":
                prediction_nfs["mu_Wminus_tot"] = sum_minus / priors["Wminus"]
                prediction_nfs["mu_Rc"] = sum_plus / sum_minus
                for i in range(1, 5):
                    prediction_nfs[f"mu_Wminus_rel_{i}"] = (gr_theory["minus"].GetY()[i - 1] / sum_minus) / (priors[f"Wminus_{i}"] / priors["Wminus"])
                    prediction_nfs[f"mu_Wplus_rel_{i}"] = (gr_theory["plus"].GetY()[i - 1] / sum_plus) / (priors[f"Wplus_{i}"] / priors["Wplus"])
            elif options.basis == "absolute":
                for i in range(1, 6):
                    prediction_nfs[f"mu_Wminus_{i}"] = gr_theory["minus"].GetY()[i - 1] / priors[f"Wminus_{i}"]
                    prediction_nfs[f"mu_Wplus_{i}"] = gr_theory["plus"].GetY()[i - 1] / priors[f"Wplus_{i}"]

            # print theory errors
            jobName = f"WCharm_lep_obs_OSSS_complete_{suffix}{var}_{prediction}"
            quickFit += " -p "
            print(f"\n============ normalization factors ({prediction}) ============")
            for x in NFS:
                print(f"{x:15s} {prediction_nfs[x]:.4f}")
                quickFit += f"{x}={prediction_nfs[x]},"
            quickFit += f"mu_Top -o {jobName}.root"
            outfile.write(f"{{\"job_name\": \"{jobName}\", \"cmd\": \"{quickFit}\"}}\n")

            # theory systematics
            if options.sys:
                # PDF theory prediction error UP (assume all correlated)
                sum_minus_up = sum([(gr_theory['minus'].GetY()[i] + gr_theory['minus'].GetEYhigh()[i]) for i in range(gr_theory['minus'].GetN())])
                sum_plus_up = sum([(gr_theory['plus'].GetY()[i] + gr_theory['plus'].GetEYhigh()[i]) for i in range(gr_theory['plus'].GetN())])
                nfs_up = {}
                nfs_up["mu_Wminus_tot"] = sum_minus_up / priors["Wminus"]
                nfs_up["mu_Rc"] = sum_plus_up / sum_minus_up
                for i in range(1, 5):
                    nfs_up[f"mu_Wminus_rel_{i}"] = prediction_nfs[f"mu_Wminus_rel_{i}"] * (1 + gr_theory_rel["minus"].GetEYhigh()[i - 1] / gr_theory_rel["minus"].GetY()[i - 1])
                    nfs_up[f"mu_Wplus_rel_{i}"] = prediction_nfs[f"mu_Wplus_rel_{i}"] * (1 + gr_theory_rel["plus"].GetEYhigh()[i - 1] / gr_theory_rel["plus"].GetY()[i - 1])
                for x in nfs_up:
                    nfs_up[x] = nfs_up[x] - prediction_nfs[x]

                # PDF theory prediction error DN (assume all correlated)
                sum_minus_dn = sum([(gr_theory['minus'].GetY()[i] - gr_theory['minus'].GetEYlow()[i]) for i in range(gr_theory['minus'].GetN())])
                sum_plus_dn = sum([(gr_theory['plus'].GetY()[i] - gr_theory['plus'].GetEYlow()[i]) for i in range(gr_theory['plus'].GetN())])
                nfs_dn = {}
                nfs_dn["mu_Wminus_tot"] = sum_minus_dn / priors["Wminus"]
                nfs_dn["mu_Rc"] = sum_plus_dn / sum_minus_dn
                for i in range(1, 5):
                    nfs_dn[f"mu_Wminus_rel_{i}"] = prediction_nfs[f"mu_Wminus_rel_{i}"] * (1 + gr_theory_rel["minus"].GetEYlow()[i - 1] / gr_theory_rel["minus"].GetY()[i - 1])
                    nfs_dn[f"mu_Wplus_rel_{i}"] = prediction_nfs[f"mu_Wplus_rel_{i}"] * (1 + gr_theory_rel["plus"].GetEYlow()[i - 1] / gr_theory_rel["plus"].GetY()[i - 1])
                for x in nfs_dn:
                    nfs_dn[x] = prediction_nfs[x] - nfs_dn[x]

                # QCD scale error file
                gr_qcd = {}
                gr_qcd_rel = {}
                for lep, meson_charge in zip(["minus", "plus"], ["plus", "minus"]):
                    if options.decay == "Dplus":
                        f_qcd = ROOT.TFile(os.path.join(PREDICTIONS_DIR, f"TheoryScaleUncert_{name}_W{lep}D{meson_charge}.root"))
                    else:
                        f_qcd = ROOT.TFile(os.path.join(PREDICTIONS_DIR, f"TheoryScaleUncert_{name}_W{lep}Dstar{meson_charge}.root"))
                    gr_qcd[lep] = f_qcd.Get(f"mu_{lep}_{options.decay}_{name}__fractionalErr_cross")
                    gr_qcd_rel[lep] = f_qcd.Get(f"mu_{lep}_{options.decay}_{name}__fractionalErr_norm")
                    assert gr_qcd[lep], (f"mu_{lep}_{options.decay}_{name}__fractionalErr_cross", prediction)

                # QCD theory prediction error UP (assume all correlated)
                sum_minus_up = sum([(gr_theory['minus'].GetY()[i] * (1 + gr_qcd['minus'].GetEYhigh()[i])) for i in range(gr_theory['minus'].GetN())])
                sum_plus_up = sum([(gr_theory['plus'].GetY()[i] * (1 + gr_qcd['plus'].GetEYhigh()[i])) for i in range(gr_theory['plus'].GetN())])
                nfs_qcd_up = {}
                nfs_qcd_up["mu_Wminus_tot"] = sum_minus_up / priors["Wminus"]
                nfs_qcd_up["mu_Rc"] = sum_plus_up / sum_minus_up
                for i in range(1, 5):
                    nfs_qcd_up[f"mu_Wminus_rel_{i}"] = prediction_nfs[f"mu_Wminus_rel_{i}"] * (1 + gr_qcd_rel["minus"].GetEYhigh()[i - 1])
                    nfs_qcd_up[f"mu_Wplus_rel_{i}"] = prediction_nfs[f"mu_Wplus_rel_{i}"] * (1 + gr_qcd_rel["plus"].GetEYhigh()[i - 1])
                for x in nfs_qcd_up:
                    nfs_qcd_up[x] = nfs_qcd_up[x] - prediction_nfs[x]

                # QCD theory prediction error DN (assume all correlated)
                sum_minus_dn = sum([(gr_theory['minus'].GetY()[i] * (1 - gr_qcd['minus'].GetEYlow()[i])) for i in range(gr_theory['minus'].GetN())])
                sum_plus_dn = sum([(gr_theory['plus'].GetY()[i] * (1 - gr_qcd['plus'].GetEYlow()[i])) for i in range(gr_theory['plus'].GetN())])
                nfs_qcd_dn = {}
                nfs_qcd_dn["mu_Wminus_tot"] = sum_minus_dn / priors["Wminus"]
                nfs_qcd_dn["mu_Rc"] = sum_plus_dn / sum_minus_dn
                for i in range(1, 5):
                    nfs_qcd_dn[f"mu_Wminus_rel_{i}"] = prediction_nfs[f"mu_Wminus_rel_{i}"] * (1 - gr_qcd_rel["minus"].GetEYlow()[i - 1])
                    nfs_qcd_dn[f"mu_Wplus_rel_{i}"] = prediction_nfs[f"mu_Wplus_rel_{i}"] * (1 - gr_qcd_rel["plus"].GetEYlow()[i - 1])
                for x in nfs_qcd_dn:
                    nfs_qcd_dn[x] = prediction_nfs[x] - nfs_qcd_dn[x]

                # covariance matrix for prediction
                # assume fully correlated and symmetrize
                cov_pred_pdf = ROOT.TMatrixD(len(NFS), len(NFS))
                for i in range(len(NFS)):
                    for j in range(len(NFS)):
                        corr = 1.0
                        if i != j:
                            corr = 0.0
                        sigma_i = (nfs_up[NFS[i]] - nfs_dn[NFS[i]]) / 2.
                        sigma_j = (nfs_up[NFS[j]] - nfs_dn[NFS[j]]) / 2.
                        cov_pred_pdf[i][j] = sigma_i * sigma_j * corr

                cov_pred_qcd = ROOT.TMatrixD(len(NFS), len(NFS))
                for i in range(len(NFS)):
                    for j in range(len(NFS)):
                        corr = 1.0
                        if i != j:
                            corr = 0.0
                        sigma_i = (nfs_qcd_up[NFS[i]] - nfs_qcd_dn[NFS[i]]) / 2.
                        sigma_j = (nfs_qcd_up[NFS[j]] - nfs_qcd_dn[NFS[j]]) / 2.
                        cov_pred_qcd[i][j] = sigma_i * sigma_j * corr

                print(f"\n============ PDF Unc. cov for {prediction} ============")
                cov_pred_pdf.Print()

                print(f"\n============ QCD Scale Unc. cov for {prediction} ============")
                cov_pred_qcd.Print()

            # invert
            cov_total = cov.Clone(f"{cov.GetName()}_{prediction}")
            if options.sys:
                cov_total += cov_pred_pdf
                cov_total += cov_pred_qcd
            icov = cov_total.Invert()

            # calculate Chi2
            chi2 = 0
            chi2_matrix = ROOT.TMatrixD(len(NFS), len(NFS))
            for i in range(len(NFS)):
                for j in range(len(NFS)):
                    chi2_i_j = (prediction_nfs[NFS[i]] - nfs[NFS[i]]) * icov[i][j] * (prediction_nfs[NFS[j]] - nfs[NFS[j]])
                    chi2_matrix[i][j] = chi2_i_j
                    chi2 += chi2_i_j
            print(f"\nChi2: {chi2}, Prob: {100 * ROOT.TMath.Prob(chi2, 10)}")
            chi2_dict[var][prediction] = chi2

            # NLL probability
            if options.quickfit_input:
                f = ROOT.TFile(os.path.join(options.input, options.quickfit_input, jobName, f"{jobName}.root"), "READ")
                t = f.Get("nllscan")
                for e in t:
                    nll = e.nll
                    break
                print(nll, nll_default)
                nll_dict[var][prediction] = 2 * (nll - nll_default)

    # print final chi2
    for var in options.vars.split(","):
        print(f"\n===== {var} =====")
        for prediction in THEORY:
            if options.quickfit_input:
                print(f"{prediction:22s}\t{100 * ROOT.TMath.Prob(chi2_dict[var][prediction], 10):.6f}\t{100 * ROOT.TMath.Prob(nll_dict[var][prediction], 10):.6f}")
            else:
                print(f"{prediction:22s}\t{100 * ROOT.TMath.Prob(chi2_dict[var][prediction], 10):.6f}")

    # close out file
    outfile.close()


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-i', '--input',
                      action="store", dest="input",
                      help="input folder with the fit results")
    parser.add_option('-s', '--sys',
                      action="store_true", dest="sys",
                      help="add theory uncertainties")
    parser.add_option('-d', '--decay',
                      action="store", dest="decay", default = "Dplus",
                      help="Decay mode (defaults to Dplus)")
    parser.add_option('-v', '--vars',
                      action="store", dest="vars",
                      help="comma sepparated list of variables",
                      default="pt,eta")
    parser.add_option('-b', '--basis',
                      action="store", dest="basis",
                      help="fit POI basis (e.g. default, absolute)",
                      default="default")
    parser.add_option('-q', '--quickfit-input',
                      action="store", dest="quickfit_input",
                      help="quickfit input folder name",
                      default="")

    # parse input arguments
    options, _ = parser.parse_args()

    # do the plotting
    main(options)
