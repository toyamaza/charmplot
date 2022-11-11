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

if not os.path.isdir("plots_sys"):
    os.makedirs("plots_sys")


def RGB(string):
    return tuple(int(string[i + 1:i + 3], 16) / 255. for i in (0, 2, 4))


colors = []
colors += [ROOT.TColor(10000, *RGB("#d3d3d3"), "10000")]
colors += [ROOT.TColor(10001, *RGB("#4FFFA1"), "10001")]
colors += [ROOT.TColor(10002, *RGB("#6A5CFF"), "10002")]
colors += [ROOT.TColor(10003, *RGB("#DFFF4F"), "10003")]
colors += [ROOT.TColor(10004, *RGB("#B32584"), "10004")]
colors += [ROOT.TColor(10005, *RGB("#FF42C1"), "10005")]
colors += [ROOT.TColor(10006, *RGB("#99B325"), "10006")]

CHANNELS = ["dplus", "dstar"]
DPLUS_FOLDER = "/global/cfs/cdirs/atlas/wcharm/TRExFitter/Output/Dplus_2022_07_26_fullRanking_v2"
DSTAR_FOLDER = "/global/cfs/cdirs/atlas/wcharm/TRExFitter/Output/Dstar_2022_08_08_fullRanking_v2"


STYLE = {
    "Luminosity": [ROOT.kGray + 3, 1, 12],
    "Signal modeling": [10005, 1, 7],
    "Muon reconstruction": [ROOT.kBlue - 9, 6, 5],
    "Electron reconstruction": [ROOT.kRed - 9, 6, 5],
    "Multijet background": [ROOT.kGreen - 8, 6, 5],
    "Finite size of MC samples": [10006, 1, 6],
    "Background modeling": [10002, 2, 4],
    "Jet and missing energy": [10001, 7, 4],
    "SV reconstruction": [10004, 1, 3],
    "Signal branching ratio": [10003, 1, 3],
}

LEGEND = {
    "SV reconstruction": "",
    "Signal modeling": "",
    "Signal branching ratio": "",
    "Finite size of MC samples": "Finite size of MC",
    "Background modeling": "",
    "Jet and missing energy": "Jet and E_{T}^{miss}",
    "Luminosity": "",
    "Muon reconstruction": "",
    "Electron reconstruction": "",
    "Multijet background": "",
}


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
        return "Signal modeling"
    else:
        return "Background modeling"


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
            sorted_groups = sorted(size.items(), key=lambda item: item[1], reverse=True)
            for g, _ in sorted_groups:
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

            # make plots for Fede
            for charge in ["minus", "plus"]:

                # fill histograms
                h_tot = ROOT.TH1D(f"h_tot_{c}_{charge}_{var}", f"h_tot_{c}_{charge}_{var}", 5, 0, 5)
                h_stat = ROOT.TH1D(f"h_stat_{c}_{charge}_{var}", f"h_stat_{c}_{charge}_{var}", 5, 0, 5)
                h_map = {}
                for g, _ in sorted_groups:
                    name = g.replace(" ", "_")
                    h_map[g] = ROOT.TH1D(f"h_{name}_{c}_{charge}_{var}", f"h_{name}_{c}_{charge}_{var}", 5, 0, 5)
                    h_map[g].SetLineColor(STYLE[g][0])
                    h_map[g].SetLineStyle(STYLE[g][1])
                    h_map[g].SetLineWidth(STYLE[g][2])

                for i in range(1, 6):
                    if charge == "minus":
                        par = PARS_DIFF[c][i - 1]
                        par_stat = PARS_STAT_DIFF[c][i - 1]
                    else:
                        par = PARS_DIFF[c][i - 1 + 5]
                        par_stat = PARS_STAT_DIFF[c][i - 1 + 5]
                    h_tot.SetBinContent(i, 100 * 0.5 * (par.getErrorHi() - par.getErrorLo()) / par.getVal())
                    h_stat.SetBinContent(i, 100 * 0.5 * (par_stat.getErrorHi() - par_stat.getErrorLo()) / par_stat.getVal())
                    POI = f"mu_W{charge}_{i}"
                    for g, _ in sorted_groups:
                        if g == "Luminosity":
                            h_map[g].SetBinContent(i, 1.7)
                        else:
                            err = 100 * (ERRORS[c][g][POI][0] + ERRORS[c][g][POI][1]) / 2.
                            if err > 0.1:
                                h_map[g].SetBinContent(i, err)
                            else:
                                h_map[g].SetBinContent(i, 1e-6)

                # histogram styles
                if var == "pt":
                    h_tot.GetXaxis().SetTitle("#it{p}_{T}^{D} [GeV]")
                elif var == "eta":
                    h_tot.GetXaxis().SetTitle("|#eta(#it{l})|")
                h_tot.GetYaxis().SetTitle("Relative Uncertaitny [%]")
                h_tot.SetMinimum(0.1)
                h_tot.SetMaximum(99)
                h_tot.GetYaxis().SetNoExponent()
                h_tot.SetFillColor(10000)
                h_tot.SetLineColor(ROOT.kGray + 1)
                h_tot.GetXaxis().SetTitleOffset(h_tot.GetXaxis().GetTitleOffset() * 1.1)
                h_tot.GetXaxis().SetLabelSize(1.5 * h_tot.GetXaxis().GetLabelSize())
                h_tot.GetYaxis().SetLabelOffset(h_tot.GetLabelOffset() * 1.2)
                if var == "pt":
                    h_tot.GetXaxis().SetLabelOffset(2.5 * h_tot.GetXaxis().GetLabelOffset())
                elif var == "eta":
                    h_tot.GetXaxis().SetLabelOffset(3.0 * h_tot.GetXaxis().GetLabelOffset())
                if var == "pt":
                    h_tot.GetXaxis().SetBinLabel(1, "8-12")
                    h_tot.GetXaxis().SetBinLabel(2, "12-20")
                    h_tot.GetXaxis().SetBinLabel(3, "20-40")
                    h_tot.GetXaxis().SetBinLabel(4, "40-80")
                    h_tot.GetXaxis().SetBinLabel(5, "#geq 80")
                elif var == "eta":
                    h_tot.GetXaxis().SetBinLabel(1, "0.0-0.5")
                    h_tot.GetXaxis().SetBinLabel(2, "0.5-1.0")
                    h_tot.GetXaxis().SetBinLabel(3, "1.0-1.5")
                    h_tot.GetXaxis().SetBinLabel(4, "1.5-2.0")
                    h_tot.GetXaxis().SetBinLabel(5, "2.0-2.5")
                # h_tot.GetXaxis().LabelsOption("v")
                h_tot.GetXaxis().ChangeLabel(1, 25)
                h_tot.GetXaxis().ChangeLabel(2, 25)
                h_tot.GetXaxis().ChangeLabel(3, 25)
                h_tot.GetXaxis().ChangeLabel(4, 25)
                h_tot.GetXaxis().ChangeLabel(5, 25)
                h_stat.SetLineColor(ROOT.kBlack)
                h_stat.SetLineWidth(4)
                h_stat.SetLineStyle(9)

                # legend
                N = 1 + len(sorted_groups)
                leg = ROOT.TLegend(0.17, 0.860 - N * (0.014), 0.94, 0.860)
                leg.SetNColumns(2)
                leg.SetBorderSize(0)
                leg.SetFillColor(0)
                leg.SetFillStyle(0)
                leg.SetTextSize(24)
                leg.SetTextFont(43)

                # plot
                canv = ROOT.TCanvas(f"sys_plot_{c}_{charge}_{var}", f"sys_plot_{c}_{charge}_{var}", 800, 1000)
                canv.SetLogy()
                h_tot.Draw("hist")
                leg.AddEntry(h_tot, "Total", "f")
                leg.AddEntry(h_stat, "Statistical", "l")
                for g in LEGEND:
                    if not LEGEND[g]:
                        leg.AddEntry(h_map[g], g, "l")
                    else:
                        leg.AddEntry(h_map[g], LEGEND[g], "l")
                for g in STYLE:
                    h_map[g].Draw("same")
                h_stat.Draw("same")

                # ATLAS label
                ROOT.ATLASLabel(0.18, 0.91, "Internal", 1, 0.04)
                ROOT.myText(0.18, 0.875, 1, "#sqrt{s} = 13 TeV, 139.0 fb^{-1}", 0.04)
                if c == "dplus":
                    if charge == "minus":
                        ROOT.myText(0.60, 0.875, 1, "W^{-}+D^{+}(#rightarrowK#pi#pi)", 0.04)
                    else:
                        ROOT.myText(0.60, 0.875, 1, "W^{+}+D^{-}(#rightarrowK#pi#pi)", 0.04)
                elif c == "dstar":
                    if charge == "minus":
                        ROOT.myText(0.60, 0.875, 1, "W^{-}+D^{*+}(#rightarrow(K#pi)#pi)", 0.04)
                    else:
                        ROOT.myText(0.60, 0.875, 1, "W^{+}+D^{*-}(#rightarrow(K#pi)#pi)", 0.04)
                leg.Draw()

                ROOT.gPad.RedrawAxis()
                canv.Print(f"plots_sys/{c}_{charge}_{var}.pdf")


if __name__ == "__main__":

    # run
    main()
