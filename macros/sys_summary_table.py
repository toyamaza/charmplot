#!/usr/bin/env python
import os
import ROOT
import yaml

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasLabels.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasUtils.C"))
ROOT.SetAtlasStyle()

CHANNELS = ["dplus", "dstar"]
DPLUS_FOLDER = "/global/cfs/cdirs/atlas/wcharm/TRExFitter/Output/Dplus_2022_07_26_fullRanking_v2"
DSTAR_FOLDER = "/global/cfs/cdirs/atlas/wcharm/TRExFitter/Output/Dstar_2022_08_08_fullRanking_v2"


def get_sys_group(name):
    if name.startswith("gamma_") or name.startswith("Stat") or name.startswith("MC_Stat"):
        return "Finite size of MC samples"
    elif name.startswith("MM_"):
        return "Multijet background"
    elif name.startswith("MUON_"):
        return "Muon reconstruction"
    elif name.startswith("Lumi"):
        return "Luminosity"
    elif name.startswith("BranchingRatio"):
        return "Signal branching ratio"
    elif name.startswith("EL_") or name.startswith("EG_"):
        return "Electron reconstruction"
    elif name.startswith("JET_") or name.startswith("FT_") or name.startswith("MET_"):
        return "Jet and missing energy"
    elif name.startswith("TRK_EFF_") or name.startswith("DMESON_MASS_") or name.startswith("DMESON_RESO_"):
        return "SV reconstruction"
    elif name.startswith("FID_EFF_") or name.startswith("Matched_Sherpa_MG_2P"):
        return "Signal modelling"
    else:
        return "Background modelling"


def main():

    for var in ["pt", "eta"]:

        ERRORS = {}
        ERRORS_DIFF = {}
        ERRORS_DIFF_NORM = {}
        PARS = {}
        PARS_STAT = {}
        PARS_DIFF = {c: [] for c in CHANNELS}
        PARS_STAT_DIFF = {c: [] for c in CHANNELS}
        PARS_DIFF_NORM = {c: [] for c in CHANNELS}
        PARS_STAT_DIFF_NORM = {c: [] for c in CHANNELS}
        size = {}

        for channel in CHANNELS:

            errors = {}
            errors_diff = {}
            errors_diff_norm = {}
            ERRORS[channel] = errors
            ERRORS_DIFF[channel] = errors_diff
            ERRORS_DIFF_NORM[channel] = errors_diff_norm

            # folder
            if channel == "dplus":
                FOLDER = DPLUS_FOLDER
            elif channel == "dstar":
                FOLDER = DSTAR_FOLDER

            # fit results
            obs_fit = f"WCharm_lep_obs_OSSS_complete_{var}"
            obs_fit2 = f"WCharm_lep_obs_OSSS_complete2_{var}"
            obs_fit3 = f"WCharm_lep_obs_OSSS_complete_alt_{var}"
            f_result = ROOT.TFile(os.path.join(FOLDER, obs_fit, "Fits", f"{obs_fit}.root"), "READ")
            f_resul2 = ROOT.TFile(os.path.join(FOLDER, obs_fit2, "Fits", f"{obs_fit2}.root"), "READ")
            f_resul3 = ROOT.TFile(os.path.join(FOLDER, obs_fit3, "Fits", f"{obs_fit3}.root"), "READ")
            fr = f_result.Get("nll_simPdf_newasimovData_with_constr")
            fr2 = f_resul2.Get("nll_simPdf_newasimovData_with_constr")
            fr3 = f_resul3.Get("nll_simPdf_newasimovData_with_constr")

            # fit results (stat. only)
            f_result_stat = ROOT.TFile(os.path.join(FOLDER, obs_fit, "Fits", f"{obs_fit}_statOnly.root"), "READ")
            f_resul2_stat = ROOT.TFile(os.path.join(FOLDER, obs_fit2, "Fits", f"{obs_fit2}_statOnly.root"), "READ")
            f_resul3_stat = ROOT.TFile(os.path.join(FOLDER, obs_fit3, "Fits", f"{obs_fit3}_statOnly.root"), "READ")
            fr_stat = f_result_stat.Get("nll_simPdf_newasimovData_with_constr")
            fr2_stat = f_resul2_stat.Get("nll_simPdf_newasimovData_with_constr")
            fr3_stat = f_resul3_stat.Get("nll_simPdf_newasimovData_with_constr")

            # print
            par_Wminus = fr.floatParsFinal().find("mu_Wminus_tot")
            par_Wplus = fr2.floatParsFinal().find("mu_Wplus_tot")
            par_Rc = fr.floatParsFinal().find("mu_Rc")
            par_Wminus_stat = fr_stat.floatParsFinal().find("mu_Wminus_tot")
            par_Wplus_stat = fr2_stat.floatParsFinal().find("mu_Wplus_tot")
            par_Rc_stat = fr_stat.floatParsFinal().find("mu_Rc")

            # inclusive and Rc
            PARS[channel] = [par_Wminus, par_Wplus, par_Rc]
            PARS_STAT[channel] = [par_Wminus_stat, par_Wplus_stat, par_Rc_stat]

            # differential bins
            for charge in ["minus", "plus"]:
                for i in range(1, 6):
                    PARS_DIFF[channel] += [fr3.floatParsFinal().find(f"mu_W{charge}_{i}")]
                    PARS_STAT_DIFF[channel] += [fr3_stat.floatParsFinal().find(f"mu_W{charge}_{i}")]
                    if i < 5:
                        PARS_DIFF_NORM[channel] += [fr.floatParsFinal().find(f"mu_W{charge}_rel_{i}")]
                        PARS_STAT_DIFF_NORM[channel] += [fr_stat.floatParsFinal().find(f"mu_W{charge}_rel_{i}")]
                    else:
                        PARS_DIFF_NORM[channel] += [fr2.floatParsFinal().find(f"mu_W{charge}_rel_{i}")]
                        PARS_STAT_DIFF_NORM[channel] += [fr2_stat.floatParsFinal().find(f"mu_W{charge}_rel_{i}")]

            # print inclusive and Rc
            for par, par_stat in zip(PARS[channel], PARS_STAT[channel]):
                par.Print()
                par_stat.Print()

            # print differential
            for par, par_stat in zip(PARS_DIFF[channel], PARS_STAT_DIFF[channel]):
                par.Print()
                par_stat.Print()

            for par, par_stat in zip(PARS_DIFF_NORM[channel], PARS_STAT_DIFF_NORM[channel]):
                par.Print()
                par_stat.Print()

            # get ranking groups
            RANKING = {}
            groups = {}

            # loop
            for POI in ["mu_Wminus_tot", "mu_Rc", "mu_Wplus_tot"]:
                fit = obs_fit
                if POI == "mu_Wplus_tot":
                    fit = obs_fit2
                with open(os.path.join(FOLDER, fit, f"Ranking_{POI}.yaml"), 'r') as stream:
                    ranking = yaml.safe_load(stream)
                    RANKING[POI] = ranking
                    for NP in ranking:
                        name = NP["Name"]
                        group = get_sys_group(name)
                        if group not in groups:
                            groups[group] = set()
                        groups[group].add(name)

            for POI in [f"mu_Wminus_{i}" for i in range(1, 6)] + [f"mu_Wplus_{i}" for i in range(1, 6)]:
                fit = obs_fit3
                with open(os.path.join(FOLDER, fit, f"Ranking_{POI}.yaml"), 'r') as stream:
                    ranking = yaml.safe_load(stream)
                    RANKING[POI] = ranking
                    for NP in ranking:
                        name = NP["Name"]
                        group = get_sys_group(name)
                        if group not in groups:
                            groups[group] = set()
                        groups[group].add(name)

            for POI in [f"mu_Wminus_rel_{i}" for i in range(1, 6)] + [f"mu_Wplus_rel_{i}" for i in range(1, 6)]:
                fit = obs_fit
                if "rel_5" in POI:
                    fit = obs_fit2
                with open(os.path.join(FOLDER, fit, f"Ranking_{POI}.yaml"), 'r') as stream:
                    ranking = yaml.safe_load(stream)
                    RANKING[POI] = ranking
                    for NP in ranking:
                        name = NP["Name"]
                        group = get_sys_group(name)
                        if group not in groups:
                            groups[group] = set()
                        groups[group].add(name)

            for k, v in groups.items():
                print("\n=======================")
                print(f"group {k} has {len(v)} systematics")
                for sys in sorted(groups[k]):
                    print(sys)

            # sum up rankings in groups
            all_POIs = ["mu_Wminus_tot", "mu_Wplus_tot", "mu_Rc"]
            all_POIs += [f"mu_Wminus_{i}" for i in range(1, 6)] + [f"mu_Wplus_{i}" for i in range(1, 6)]
            all_POIs += [f"mu_Wminus_rel_{i}" for i in range(1, 6)] + [f"mu_Wplus_rel_{i}" for i in range(1, 6)]
            all_pars = PARS[channel] + PARS_DIFF[channel] + PARS_DIFF_NORM[channel]
            print("\n=======================")
            for k, v in groups.items():
                errors[k] = {}
                if k not in size:
                    size[k] = 0
                for POI, par in zip(all_POIs, all_pars):
                    err_up = 0
                    err_dn = 0
                    ranking = RANKING[POI]
                    for sys in v:
                        NP_impact = [x for x in ranking if x['Name'] == sys]
                        if len(NP_impact):
                            NP = NP_impact[0]
                            for key in ["POIup", "POIdown"]:
                                if NP[key] > 0:
                                    err_up += NP[key]**2
                                else:
                                    err_dn += NP[key]**2
                        else:
                            print(f"WARNING: NP {sys} not found for POI {POI}")
                    print(f"{f'Group {k}, POI {POI}':40s} +{err_up**0.5 / par.getVal():.3f} -{err_dn**0.5 / par.getVal():.3f}")
                    errors[k][POI] = [err_up**0.5 / par.getVal(), err_dn**0.5 / par.getVal()]
                    if POI in ["mu_Wminus_tot", "mu_Wplus_tot", "mu_Rc"]:
                        size[k] += err_up + err_dn

        # print LaTeX format (inclusive)
        print(f"\n============== {var} ==============")
        for g, _ in sorted(size.items(), key=lambda item: item[1], reverse=True):
            errs = []
            for POI in ["mu_Wminus_tot", "mu_Wplus_tot", "mu_Rc"]:
                errs += [(ERRORS["dplus"][g][POI][0] + ERRORS["dplus"][g][POI][1]) / 2.]
            for POI in ["mu_Wminus_tot", "mu_Wplus_tot", "mu_Rc"]:
                errs += [(ERRORS["dstar"][g][POI][0] + ERRORS["dstar"][g][POI][1]) / 2.]
            print(f"    {g:30s}{' '.join([f'& {100 * err:.1f}' for err in errs])} \\\\")
        errs = [par.getErrorHi() / par.getVal() for par in PARS_STAT["dplus"]]
        errs += [par.getErrorHi() / par.getVal() for par in PARS_STAT["dstar"]]
        print(f"    {'Data statistical uncertainty':30s}{' '.join([f'& {100 * err:.1f}' for err in errs])} \\\\")
        print(r"    \midrule")
        errs = [0.5 * (par.getErrorHi() - par.getErrorLo()) / par.getVal() for par in PARS["dplus"]]
        errs += [0.5 * (par.getErrorHi() - par.getErrorLo()) / par.getVal() for par in PARS["dstar"]]
        print(f"    {'Total':30s}{' '.join([f'& {100 * err:.1f}' for err in errs])} \\\\")

        # print LaTeX format (differential)
        for c in CHANNELS:
            print(f"\n============= {c} {var} ===============")
            for g, _ in sorted(size.items(), key=lambda item: item[1], reverse=True):
                errs = []
                errs_norm = []
                for POI in [f"mu_Wminus_{i}" for i in range(1, 6)] + [f"mu_Wplus_{i}" for i in range(1, 6)]:
                    errs += [(ERRORS[c][g][POI][0] + ERRORS[c][g][POI][1]) / 2.]
                for POI in [f"mu_Wminus_rel_{i}" for i in range(1, 6)] + [f"mu_Wplus_rel_{i}" for i in range(1, 6)]:
                    errs_norm += [(ERRORS[c][g][POI][0] + ERRORS[c][g][POI][1]) / 2.]
                print(f"    {g:30s}{' '.join([f'& {100 * err:.1f} ({100 * err2:.1f})' for (err, err2) in zip(errs, errs_norm)])} \\\\")
            errs = [par.getErrorHi() / par.getVal() for par in PARS_STAT_DIFF[c]]
            errs_norm = [par.getErrorHi() / par.getVal() for par in PARS_STAT_DIFF_NORM[c]]
            print(f"    {'Data statistical uncertainty':30s}{' '.join([f'& {100 * err:.1f} ({100 * err2:.1f})' for (err, err2) in zip(errs, errs_norm)])} \\\\")
            print(r"    \midrule")
            errs = [0.5 * (par.getErrorHi() - par.getErrorLo()) / par.getVal() for par in PARS_DIFF[c]]
            errs_norm = [0.5 * (par.getErrorHi() - par.getErrorLo()) / par.getVal() for par in PARS_DIFF_NORM[c]]
            print(f"    {'Total':30s}{' '.join([f'& {100 * err:.1f} ({100 * err2:.1f})' for (err, err2) in zip(errs, errs_norm)])} \\\\")


if __name__ == "__main__":

    # run
    main()
