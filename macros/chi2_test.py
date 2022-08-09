#!/usr/bin/env python
import os
import ROOT
import yaml

PREDICTIONS_DIR = "/global/cfs/cdirs/atlas/wcharm/Rivet/v1/processed"
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


def decode_rel_index(i):
    if i < 4:
        index_i = i
    elif i == 4:
        index_i = 11
    elif i < 9:
        index_i = i
    elif i == 9:
        index_i = 12
    return index_i

def main(options):

    # save chi2
    chi2_dict = {}
    nll_dict = {}

    # quickFit command list
    outfile = open("task_list_NLL_test.txt", "w")

    # loop over stuff
    for var in options.vars.split(","):

        # normalization factors basis
        if options.basis == "relative":
            #       0 -> 0            1 -> 1            2 -> 2            3 -> 3            4 -> 11
            NFS = ["mu_Wplus_rel_1", "mu_Wplus_rel_2", "mu_Wplus_rel_3", "mu_Wplus_rel_4", "mu_Wminus_tot",
            #       5 -> 5 (0 + 5)     6 -> 6 (1 + 5)     7 -> 7             8 -> 8             9 -> 12
                   "mu_Wminus_rel_1", "mu_Wminus_rel_2", "mu_Wminus_rel_3", "mu_Wminus_rel_4", "mu_Rc"]
            suffix = ""
        elif options.basis == "absolute":
            NFS = ["mu_Wplus_1", "mu_Wplus_2", "mu_Wplus_3", "mu_Wplus_4", "mu_Wplus_5",
                   "mu_Wminus_1", "mu_Wminus_2", "mu_Wminus_3", "mu_Wminus_4", "mu_Wminus_5"]
            suffix = "alt_"

        # cross section priors
        f_prior = ROOT.TFile(os.path.join(PRIORS_DIR, f"fid_eff_{var}_{options.decay.lower()}", "unfolding.root"), "READ")
        h_minus = f_prior.Get(f"Sherpa2211_WplusD_OS-SS_lep_minus_{options.decay}_Kpipi_truth_differential_{var}")
        h_plus = f_prior.Get(f"Sherpa2211_WplusD_OS-SS_lep_plus_{options.decay}_Kpipi_truth_differential_{var}")

        # Multiply scale factor by 0.677 for D*
        br = 1.0
        if options.decay == "Dstar":
            br = 0.677

        sf = 1.0 * br
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

        # covariance matrix from the ranking plots
        cov_ranking = ROOT.TMatrixD(len(NFS), len(NFS))
        if os.path.isdir(os.path.join(options.input, f"WCharm_lep_obs_OSSS_complete_{suffix}{var}", "Covariances")):
            with open(os.path.join(options.input, f"WCharm_lep_obs_OSSS_complete_{suffix}{var}", "Covariances", "Combined_Covariance_postfit.yaml"), 'r') as stream:
                cov_dict = yaml.safe_load(stream)
                for i in range(len(NFS)):
                    for j in range(len(NFS)):
                        index_i = cov_dict[0]["parameters"].index(NFS[i])
                        index_j = cov_dict[0]["parameters"].index(NFS[j])
                        cov_ranking[i][j] = cov_dict[1]["covariance_rows"][index_i][index_j]

        print(f"\n============ Experimental cov for {var} ============")
        cov_ranking.Print()
        if options.ranking_cov:
            cov = cov_ranking
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
            if "NNPDF" in prediction:
                prediction += "_hessian"
            name = "lep_abs_eta" if var == "eta" else "D_pt"
            f_theory = ROOT.TFile(os.path.join(PREDICTIONS_DIR, "predictions.root"))
            assert f_theory

            # quickFit template
            if options.decay == "Dplus":
                quickFit = f"quickFit --minStrat 1 --minTolerance 1e-5 --savefitresult 0 --hesse 0 --saveNP 0 -f {options.input}/WCharm_lep_obs_OSSS_complete_{suffix}{var}/RooStats/WCharm_lep_obs_OSSS_complete_{suffix}{var}_combined_WCharm_lep_obs_OSSS_complete_{suffix}{var}_model.root -w combined -d obsData "
            else:
                quickFit = f"quickFit --minStrat 1 --minTolerance 1e-5 --savefitresult 0 --hesse 0 --saveNP 0 -f {options.input}/WCharm_lep_obs_OSSS_complete_{suffix}{var}/RooStats/WCharm_lep_obs_OSSS_complete{suffix}_combined_WCharm_lep_obs_OSSS_complete{suffix}_model.root -w combined -d obsData "

            # read theory histograms
            # absolute basis
            hName_base = f"{prediction}_abs_{options.decay}_{name}"
            h_theory = f_theory.Get(f"{hName_base}_pdf_central")
            h_theory_pdf_up = f_theory.Get(f"{hName_base}_pdf_up")
            h_theory_pdf_dn = f_theory.Get(f"{hName_base}_pdf_dn")
            assert (h_theory and h_theory_pdf_up and h_theory_pdf_dn), hName_base

            # relative basis
            hName_rel_base = f"{prediction}_rel_{options.decay}_{name}"
            h_theory_rel = f_theory.Get(f"{hName_rel_base}_pdf_central")
            h_theory_rel_pdf_up = f_theory.Get(f"{hName_rel_base}_pdf_up")
            h_theory_rel_pdf_dn = f_theory.Get(f"{hName_rel_base}_pdf_dn")
            assert (h_theory_rel and h_theory_rel_pdf_up and h_theory_rel_pdf_dn), hName_rel_base

            # normalization factor values for the prediction
            prediction_nfs = {}
            if options.basis == "absolute":
                for i in range(1, 6):
                    # absolute cross section
                    # POI order: 5 x W+D-, 5 x W-D+
                    prediction_nfs[f"mu_Wplus_{i}"] = h_theory.GetBinContent(i) / priors[f"Wplus_{i}"]
                    prediction_nfs[f"mu_Wminus_{i}"] = h_theory.GetBinContent(i + 5) / priors[f"Wminus_{i}"]
            elif options.basis == "relative":
                for i in range(1, 5):
                    # relative cross section
                    # POI order: 5 x W+D-, 5 x W-D+, W+ tot, W- tot, Rc
                    prediction_nfs[f"mu_Wplus_rel_{i}"] = h_theory_rel.GetBinContent(i) / (priors[f"Wplus_{i}"] / priors[f"Wplus"])
                    prediction_nfs[f"mu_Wminus_rel_{i}"] = h_theory_rel.GetBinContent(i + 5) / (priors[f"Wminus_{i}"] / priors[f"Wminus"])
                sum_minus = sum([h_theory.GetBinContent(i + 5) for i in range(1, 6)])
                sum_plus = sum([h_theory.GetBinContent(i) for i in range(1, 6)])
                prediction_nfs[f"mu_Wminus_tot"] = sum_minus / priors[f"Wminus"]
                prediction_nfs[f"mu_Rc"] = sum_plus / sum_minus

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
                # PDF theory prediction error
                nfs_up = {}
                nfs_dn = {}
                if options.basis == "absolute":
                    for i in range(1, 6):
                        nfs_up[f"mu_Wplus_{i}"] = h_theory_pdf_up.GetBinContent(i) / priors[f"Wplus_{i}"]
                        nfs_up[f"mu_Wminus_{i}"] = h_theory_pdf_up.GetBinContent(i + 5) / priors[f"Wminus_{i}"]
                        nfs_dn[f"mu_Wplus_{i}"] = h_theory_pdf_dn.GetBinContent(i) / priors[f"Wplus_{i}"]
                        nfs_dn[f"mu_Wminus_{i}"] = h_theory_pdf_dn.GetBinContent(i + 5) / priors[f"Wminus_{i}"]
                elif options.basis == "relative":
                    for i in range(1, 5):
                        nfs_up[f"mu_Wplus_rel_{i}"] = h_theory_rel_pdf_up.GetBinContent(i) / (priors[f"Wplus_{i}"] / priors[f"Wplus"])
                        nfs_up[f"mu_Wminus_rel_{i}"] = h_theory_rel_pdf_up.GetBinContent(i + 5) / (priors[f"Wminus_{i}"] / priors[f"Wminus"])
                        nfs_dn[f"mu_Wplus_rel_{i}"] = h_theory_rel_pdf_dn.GetBinContent(i) / (priors[f"Wplus_{i}"] / priors[f"Wplus"])
                        nfs_dn[f"mu_Wminus_rel_{i}"] = h_theory_rel_pdf_dn.GetBinContent(i + 5) / (priors[f"Wminus_{i}"] / priors[f"Wminus"])
                    sum_plus_up = sum([h_theory_pdf_up.GetBinContent(i) for i in range(1, 6)])
                    sum_minus_up = sum([h_theory_pdf_up.GetBinContent(i + 5) for i in range(1, 6)])
                    sum_plus_dn = sum([h_theory_pdf_dn.GetBinContent(i) for i in range(1, 6)])
                    sum_minus_dn = sum([h_theory_pdf_dn.GetBinContent(i + 5) for i in range(1, 6)])
                    nfs_up[f"mu_Wminus_tot"] = sum_minus_up / priors[f"Wminus"]
                    nfs_up[f"mu_Rc"] = sum_plus_up / sum_minus_up
                    nfs_dn[f"mu_Wminus_tot"] = sum_minus_dn / priors[f"Wminus"]
                    nfs_dn[f"mu_Rc"] = sum_plus_dn / sum_minus_dn

                # # QCD scale error file
                # h_qcd_up = {}
                # h_qcd_dn = {}
                # h_qcd_rel_up = {}
                # h_qcd_rel_dn = {}
                # if options.decay == "Dplus":
                #     f_qcd = ROOT.TFile(os.path.join(PREDICTIONS_DIR, f"ScalePredictions_{name}_WD.root"))
                # else:
                #     f_qcd = ROOT.TFile(os.path.join(PREDICTIONS_DIR, f"ScalePredictions_{name}_WDstar.root"))
                # for lep, _ in zip(["minus", "plus"], ["plus", "minus"]):
                #     h_qcd_up[lep] = f_qcd.Get(f"mu_{lep}_{options.decay}_{name}_ScaleupRatio_cross")
                #     h_qcd_dn[lep] = f_qcd.Get(f"mu_{lep}_{options.decay}_{name}_ScaledownRatio_cross")
                #     h_qcd_rel_up[lep] = f_qcd.Get(f"mu_{lep}_{options.decay}_{name}_norm_ScaleupRatio_cross")
                #     h_qcd_rel_dn[lep] = f_qcd.Get(f"mu_{lep}_{options.decay}_{name}_norm_ScaledownRatio_cross")
                #     assert h_qcd_up[lep], f"mu_{lep}_{options.decay}_{name}_ScaleupRatio_cross"

                # # QCD theory prediction error
                # nfs_qcd_up = {}
                # nfs_qcd_dn = {}
                # if options.basis == "default":
                #     sum_minus_up = sum([(h_theory['minus'].GetBinContent(i) * h_qcd_up['minus'].GetBinContent(i)) for i in range(1, h_theory['minus'].GetBinsX() + 1)])
                #     sum_minus_dn = sum([(h_theory['minus'].GetBinContent(i) * h_qcd_dn['minus'].GetBinContent(i)) for i in range(1, h_theory['minus'].GetBinsX() + 1)])
                #     sum_plus_up = sum([(h_theory['plus'].GetBinContent(i) * h_qcd_up['plus'].GetBinContent(i)) for i in range(1, h_theory['plus'].GetBinsX() + 1)])
                #     sum_plus_dn = sum([(h_theory['plus'].GetBinContent(i) * h_qcd_dn['plus'].GetBinContent(i)) for i in range(1, h_theory['plus'].GetBinsX() + 1)])
                #     nfs_qcd_up["mu_Wminus_tot"] = sum_minus_up / priors["Wminus"]
                #     nfs_qcd_dn["mu_Wminus_tot"] = sum_minus_dn / priors["Wminus"]
                #     nfs_qcd_up["mu_Rc"] = sum_plus_up / sum_minus_up
                #     nfs_qcd_dn["mu_Rc"] = sum_plus_dn / sum_minus_dn
                #     for i in range(1, 5):
                #         nfs_qcd_up[f"mu_Wminus_rel_{i}"] = prediction_nfs[f"mu_Wminus_rel_{i}"] * h_qcd_rel_up['minus'].GetBinContent(i)
                #         nfs_qcd_dn[f"mu_Wminus_rel_{i}"] = prediction_nfs[f"mu_Wminus_rel_{i}"] * h_qcd_rel_dn['minus'].GetBinContent(i)
                #         nfs_qcd_up[f"mu_Wplus_rel_{i}"] = prediction_nfs[f"mu_Wplus_rel_{i}"] * h_qcd_rel_up['plus'].GetBinContent(i)
                #         nfs_qcd_dn[f"mu_Wplus_rel_{i}"] = prediction_nfs[f"mu_Wplus_rel_{i}"] * h_qcd_rel_dn['plus'].GetBinContent(i)
                # elif options.basis == "absolute":
                #     for i in range(1, 6):
                #         nfs_qcd_up[f"mu_Wminus_{i}"] = prediction_nfs[f"mu_Wminus_{i}"] * h_qcd_up['minus'].GetBinContent(i)
                #         nfs_qcd_dn[f"mu_Wminus_{i}"] = prediction_nfs[f"mu_Wminus_{i}"] * h_qcd_dn['minus'].GetBinContent(i)
                #         nfs_qcd_up[f"mu_Wplus_{i}"] = prediction_nfs[f"mu_Wplus_{i}"] * h_qcd_up['plus'].GetBinContent(i)
                #         nfs_qcd_dn[f"mu_Wplus_{i}"] = prediction_nfs[f"mu_Wplus_{i}"] * h_qcd_dn['plus'].GetBinContent(i)

                # covariance matrix for prediction
                if options.basis == "absolute":
                    m_corr = f_theory.Get(f"{hName_base}_pdf_corr")
                    assert m_corr, f"{hName_base}_pdf_corr"
                elif options.basis == "relative":
                    m_corr = f_theory.Get(f"{hName_rel_base}_pdf_corr")
                    assert m_corr, f"{hName_rel_base}_pdf_corr"

                cov_pred_pdf = ROOT.TMatrixD(len(NFS), len(NFS))
                for i in range(len(NFS)):
                    for j in range(len(NFS)):
                        corr = 1.0
                        if i != j:
                            corr = 1.0
                        if options.basis == "absolute":
                            corr = m_corr.GetBinContent(i + 1, j + 1)
                        elif options.basis == "relative":
                            corr = m_corr.GetBinContent(decode_rel_index(i) + 1, decode_rel_index(j) + 1)
                            print("XXX ", decode_rel_index(i) + 1, decode_rel_index(j) + 1, NFS[i], NFS[j], corr)
                        sigma_i = (nfs_up[NFS[i]] - nfs_dn[NFS[i]]) / 2.
                        sigma_j = (nfs_up[NFS[j]] - nfs_dn[NFS[j]]) / 2.
                        cov_pred_pdf[i][j] = sigma_i * sigma_j * corr

                # cov_pred_qcd = ROOT.TMatrixD(len(NFS), len(NFS))
                # for i in range(len(NFS)):
                #     for j in range(len(NFS)):
                #         corr = 1.0
                #         if i != j:
                #             corr = 1.0
                #         sigma_i = (nfs_qcd_up[NFS[i]] - nfs_qcd_dn[NFS[i]]) / 2.
                #         sigma_j = (nfs_qcd_up[NFS[j]] - nfs_qcd_dn[NFS[j]]) / 2.
                #         cov_pred_qcd[i][j] = sigma_i * sigma_j * corr

                # production fraction and hadronization
                if options.basis == "absolute":
                    cov_pred_hadronization = ROOT.TMatrixD(len(NFS), len(NFS))
                    for i in range(len(NFS)):
                        for j in range(len(NFS)):
                            corr = 1.0
                            if i != j:
                                corr = 1.0
                            if options.decay == "Dplus":
                                sigma_i = (0.01**2 + 0.028**2)**(0.5)
                                sigma_j = (0.01**2 + 0.028**2)**(0.5)
                            elif options.decay == "Dstar":
                                sigma_i = (0.01**2 + 0.028**2)**(0.5)
                                sigma_j = (0.01**2 + 0.028**2)**(0.5)
                            cov_pred_hadronization[i][j] = sigma_i * sigma_j * corr

                print(f"\n============ PDF Unc. cov for {prediction} ============")
                for x in NFS:
                    print(f"{x:15s} {nfs_up[x]:.4f} {prediction_nfs[x]:.4f} {nfs_dn[x]:.4f}")
                cov_pred_pdf.Print()

                # print(f"\n============ QCD Scale Unc. cov for {prediction} ============")
                # for x in NFS:
                #     print(f"{x:15s} {nfs_qcd_up[x]:.4f} {prediction_nfs[x]:.4f} {nfs_qcd_dn[x]:.4f}")
                # cov_pred_qcd.Print()

                if options.basis == "absolute":
                    print(f"\n============ Prod Frac. cov for {prediction} ============")
                    cov_pred_hadronization.Print()

            # invert
            cov_total = cov.Clone(f"{cov.GetName()}_{prediction}")
            if options.sys:
                cov_total += cov_pred_pdf
                # cov_total += cov_pred_qcd
                if options.basis == "absolute" and options.production_fraction_sys:
                    cov_total += cov_pred_hadronization
            icov = cov_total.Invert()

            # calculate Chi2
            chi2 = 0
            chi2_matrix = ROOT.TMatrixD(len(NFS), len(NFS))
            chi2_array = []
            for i in range(len(NFS)):
                for j in range(len(NFS)):
                    print(i, j, prediction_nfs[NFS[i]] - nfs[NFS[i]], icov[i][j], prediction_nfs[NFS[j]] - nfs[NFS[j]])
                    chi2_i_j = (prediction_nfs[NFS[i]] - nfs[NFS[i]]) * icov[i][j] * (prediction_nfs[NFS[j]] - nfs[NFS[j]])
                    chi2_matrix[i][j] = chi2_i_j
                    chi2 += chi2_i_j
                    chi2_array += [(chi2_i_j, NFS[i], NFS[j])]
            chi2_dict[var][prediction.replace("_hessian", "")] = chi2
            chi2_array.sort(key=lambda x: abs(x[0]), reverse=True)
            chi2_matrix.Print()
            print(f"\nChi2: {chi2}, Prob: {100 * ROOT.TMath.Prob(chi2, 10)}\n")
            partial_chi2 = 0
            for x in chi2_array:
                partial_chi2 += x[0]
                print(f"{partial_chi2:.2f} {x[0]:.2f} {x[1]} {x[2]}")

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
    parser.add_option('-r', '--ranking-cov',
                      action="store_true", dest="ranking_cov",
                      help="use the covariance matrix from the ranking plots")
    parser.add_option('-v', '--vars',
                      action="store", dest="vars",
                      help="comma sepparated list of variables",
                      default="pt,eta")
    parser.add_option('-b', '--basis',
                      action="store", dest="basis",
                      help="fit POI basis (e.g. absolute, relative)",
                      default="absolute")
    parser.add_option('-q', '--quickfit-input',
                      action="store", dest="quickfit_input",
                      help="quickfit input folder name",
                      default="")
    parser.add_option('-p', '--production-fraction-sys',
                      action="store_true", dest="production_fraction_sys",
                      help="add production fraction systematics")

    # parse input arguments
    options, _ = parser.parse_args()

    # do the plotting
    main(options)
