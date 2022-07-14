#!/usr/bin/env python
import ROOT
import os

PREDICTIONS_DIR = "/global/cfs/cdirs/atlas/shapiro/13tev_wd_rivet/forMiha"
PRIORS_DIR = "/global/cfs/cdirs/atlas/wcharm/charmplot_output/Dmeson_2022_06_15/"
UNCERTAINTIES_DIR = "/global/cfs/cdirs/atlas/wcharm/wplusd_uncertainties"

NFS = ["mu_Rc", "mu_Wminus_rel_1", "mu_Wminus_rel_2", "mu_Wminus_rel_3", "mu_Wminus_rel_4",
       "mu_Wminus_tot", "mu_Wplus_rel_1", "mu_Wplus_rel_2", "mu_Wplus_rel_3", "mu_Wplus_rel_4"]

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

    # loop over stuff
    for var in ["eta", "pt"]:

        # cross section priors
        f_prior = ROOT.TFile(os.path.join(PRIORS_DIR, f"fid_eff_{var}_dplus", "unfolding.root"), "READ")
        h_minus = f_prior.Get(f"Sherpa2211_WplusD_OS-SS_lep_minus_Dplus_Kpipi_truth_differential_{var}")
        h_plus = f_prior.Get(f"Sherpa2211_WplusD_OS-SS_lep_plus_Dplus_Kpipi_truth_differential_{var}")

        # TODO: fix for branching ratio in W+D* forced decay sample
        br = 1.0
        priors = {
            "Wminus": h_minus.Integral() / (2. * br),
            "Wminus_1": h_minus.GetBinContent(1) / (2. * br),
            "Wminus_2": h_minus.GetBinContent(2) / (2. * br),
            "Wminus_3": h_minus.GetBinContent(3) / (2. * br),
            "Wminus_4": h_minus.GetBinContent(4) / (2. * br),
            "Wminus_5": h_minus.GetBinContent(5) / (2. * br),
            "Wplus": h_plus.Integral() / (2. * br),
            "Wplus_1": h_plus.GetBinContent(1) / (2. * br),
            "Wplus_2": h_plus.GetBinContent(2) / (2. * br),
            "Wplus_3": h_plus.GetBinContent(3) / (2. * br),
            "Wplus_4": h_plus.GetBinContent(4) / (2. * br),
            "Wplus_5": h_plus.GetBinContent(5) / (2. * br),
        }
        print(f"============ cross section priors for {var} ============")
        for key, val in priors.items():
            print(f"{key}: {val}")

        # normalization factors
        f_result = ROOT.TFile(os.path.join(options.input, f"WCharm_lep_obs_OSSS_complete_{var}", "Fits", f"WCharm_lep_obs_OSSS_complete_{var}.root"), "READ")
        fr = f_result.Get("nll_simPdf_newasimovData_with_constr")
        nfs = {nf: fr.floatParsFinal().find(nf).getVal() for nf in NFS}
        print(nfs)

        # covariance matrix
        cov = ROOT.TMatrixD(len(NFS), len(NFS))
        for i in range(len(NFS)):
            for j in range(len(NFS)):
                corr = fr.correlation(NFS[i], NFS[j])
                sigma_i = fr.floatParsFinal().find(NFS[i]).getError()
                sigma_j = fr.floatParsFinal().find(NFS[j]).getError()
                cov[i][j] = sigma_i * sigma_j * corr
        cov.Print()

        # invert
        icov = cov.Invert()
        icov.Print()

        # get prediction
        chi2_dict[var] = {}
        for prediction in THEORY:
            name = "lep_eta" if var == "eta" else "D_pt"
            lep = "minus"
            meson_charge = "plus"
            f_theory = {}
            for lep, meson_charge in zip(["minus", "plus"], ["plus", "minus"]):
                f_name = f"TheoryPredictions_{name}_W{lep}D{meson_charge}.root"
                f = ROOT.TFile(os.path.join(PREDICTIONS_DIR, f_name))
                f_theory[lep] = f
                assert f

            # read graphs
            gr_theory = {}
            for lep, meson_charge in zip(["minus", "plus"], ["plus", "minus"]):
                gr = f_theory[lep].Get(f"mu_{lep}_Dplus_{name}_{prediction}_cross")
                assert gr, f"mu_{lep}_Dplus_{name}_{prediction}_cross"
                gr_theory[lep] = gr
            sum_minus = sum([gr_theory['minus'].GetY()[i] for i in range(gr_theory['minus'].GetN())])
            sum_plus = sum([gr_theory['plus'].GetY()[i] for i in range(gr_theory['plus'].GetN())])

            # normalization factor values for the prediction
            prediction_nfs = {}
            if var == "eta":
                prediction_nfs["mu_Wminus_tot"] = (sum_minus / 2.) / priors["Wminus"]
            else:
                prediction_nfs["mu_Wminus_tot"] = sum_minus / priors["Wminus"]
            prediction_nfs["mu_Wminus_rel_1"] = (gr_theory["minus"].GetY()[0] / sum_minus) / (priors["Wminus_1"] / priors["Wminus"])
            prediction_nfs["mu_Wminus_rel_2"] = (gr_theory["minus"].GetY()[1] / sum_minus) / (priors["Wminus_2"] / priors["Wminus"])
            prediction_nfs["mu_Wminus_rel_3"] = (gr_theory["minus"].GetY()[2] / sum_minus) / (priors["Wminus_3"] / priors["Wminus"])
            prediction_nfs["mu_Wminus_rel_4"] = (gr_theory["minus"].GetY()[3] / sum_minus) / (priors["Wminus_4"] / priors["Wminus"])
            prediction_nfs["mu_Rc"] = sum_plus / sum_minus
            prediction_nfs["mu_Wplus_rel_1"] = (gr_theory["plus"].GetY()[0] / sum_plus) / (priors["Wplus_1"] / priors["Wplus"])
            prediction_nfs["mu_Wplus_rel_2"] = (gr_theory["plus"].GetY()[1] / sum_plus) / (priors["Wplus_2"] / priors["Wplus"])
            prediction_nfs["mu_Wplus_rel_3"] = (gr_theory["plus"].GetY()[2] / sum_plus) / (priors["Wplus_3"] / priors["Wplus"])
            prediction_nfs["mu_Wplus_rel_4"] = (gr_theory["plus"].GetY()[3] / sum_plus) / (priors["Wplus_4"] / priors["Wplus"])
            print(prediction_nfs)

            # calculate Chi2
            chi2 = 0
            chi2_matrix = ROOT.TMatrixD(len(NFS), len(NFS))
            for i in range(len(NFS)):
                for j in range(len(NFS)):
                    chi2_i_j = (prediction_nfs[NFS[i]] - nfs[NFS[i]]) * icov[i][j] * (prediction_nfs[NFS[j]] - nfs[NFS[j]])
                    chi2_matrix[i][j] = chi2_i_j
                    chi2 += chi2_i_j
            chi2_matrix.Print()
            print(f"Chi2: {chi2}, Prob: {ROOT.TMath.Prob(chi2, 10)}")
            chi2_dict[var][prediction] = chi2

    # print final chi2
    for var in ["pt", "eta"]:
        print(f"===== {var} =====")
        for prediction in THEORY:
            print(f"{prediction:26s}:\t{100 * ROOT.TMath.Prob(chi2_dict[var][prediction], 10):.6f}")


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-i', '--input',
                      action="store", dest="input",
                      help="input folder with the fit results")

    # parse input arguments
    options, args = parser.parse_args()

    # do the plotting
    main(options)
