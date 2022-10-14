#!/usr/bin/env python
from math import log10, floor
import os
import ROOT
import yaml

ROOT.gROOT.SetBatch(True)

DPLUS_FOLDER = "/global/cfs/cdirs/atlas/wcharm/TRExFitter/Output/Dplus_2022_08_05_v2"
DSTAR_FOLDER = "/global/cfs/cdirs/atlas/wcharm/TRExFitter/Output/Dstar_2022_08_11"

POIs_abs = [f"mu_Wplus_{i}" for i in range(1, 6)] + [f"mu_Wminus_{i}" for i in range(1, 6)]

PRUNING_THRESHOLD = 0.98

QUALIFIERS = [{'name': 'SQRT(s)', 'units': 'GeV', 'value': 13000}, {'name': 'LUMINOSITY', 'units': 'fb$^{-1}$', 'value': 139}]

CHANNELS = ["dplus", "dstar"]

# make output folder
if not os.path.isdir("hepdata"):
    os.makedirs("hepdata")


def atlas_rounding(central, up, dn=None, precision=None):
    if precision:
        central = round(central, precision)
        up = round(up, precision)
        if dn:
            dn = round(dn, precision)
        return central, up, dn
    up = round(up, -int(floor(log10(abs(up / 100.)))))
    if dn:
        dn = round(dn, -int(floor(log10(abs(dn / 100.)))))
        central = round(central, min(-int(floor(log10(abs(up / 100.)))), -int(floor(log10(abs(dn / 100.))))))
        return central, up, dn, min(-int(floor(log10(abs(up / 100.)))), -int(floor(log10(abs(dn / 100.)))))
    else:
        central = round(central, -int(floor(log10(abs(up / 100.)))))
        return central, up, -int(floor(log10(abs(up / 100.))))


def bin_edges(POI, var):
    if var == "eta":
        i = int(POI[-1])
        return (0.5 * (i - 1), 0.5 * i)


def dependable_dict(name):
    return {'header': {'name': name},
            'qualifiers': QUALIFIERS,
            'values': []}


def main():

    # main yaml file
    submission = [
        {
            "data_file": "dplus_minus.yaml",
            "description": "The 'OS-SS' W+D fiducial phase-space absolute differential cross-sections in the W-D+ channel.",
            "keywords": [
                {"name": "reactions", "values": ["P P --> W D"]},
                {"name": "observables", "values": ["SIG"]},
                {"name": "cmenergies", "values": [13000.0]},
                {"name": "phrases", "values": ["Differential Cross Section", "Cross Section", "Proton-Proton Scattering", "W Production", "D production"]},
            ],
            "location": "XXX",
            "name": "NP Impact on OS-SS W+D Cross Section (W-D+ channel)"
        },
        {
            "data_file": "dplus_plus.yaml",
            "description": "The 'OS-SS' W+D fiducial phase-space absolute differential cross-sections in the W+D- channel.",
            "keywords": [
                {"name": "reactions", "values": ["P P --> W D"]},
                {"name": "observables", "values": ["SIG"]},
                {"name": "cmenergies", "values": [13000.0]},
                {"name": "phrases", "values": ["Differential Cross Section", "Cross Section", "Proton-Proton Scattering", "W Production", "D production"]},
            ],
            "location": "XXX",
            "name": "NP Impact on OS-SS W+D Cross Section (W+D- channel)"
        },
        {
            "data_file": "dstar_minus.yaml",
            "description": "The 'OS-SS' W+D* fiducial phase-space absolute differential cross-sections in the W-D*+ channel.",
            "keywords": [
                {"name": "reactions", "values": ["P P --> W D"]},
                {"name": "observables", "values": ["SIG"]},
                {"name": "cmenergies", "values": [13000.0]},
                {"name": "phrases", "values": ["Differential Cross Section", "Cross Section", "Proton-Proton Scattering", "W Production", "D production"]},
            ],
            "location": "XXX",
            "name": "NP Impact on OS-SS W+D Cross Section (W-D*+ channel)"
        },
        {
            "data_file": "dstar_plus.yaml",
            "description": "The 'OS-SS' W+D* fiducial phase-space absolute differential cross-sections in the W+D*- channel.",
            "keywords": [
                {"name": "reactions", "values": ["P P --> W D"]},
                {"name": "observables", "values": ["SIG"]},
                {"name": "cmenergies", "values": [13000.0]},
                {"name": "phrases", "values": ["Differential Cross Section", "Cross Section", "Proton-Proton Scattering", "W Production", "D production"]},
            ],
            "location": "XXX",
            "name": "NP Impact on OS-SS W+D Cross Section (W+D*- channel)"
        },
    ]
    with open('hepdata/submission.yaml', 'w') as yaml_file:
        for table in submission:
            yaml_file.write(f"# Nuisance parameter ranking for {table['data_file'].replace('.yaml', '')}\n")
            yaml_file.write("---\n")
            yaml.dump(table, yaml_file, default_flow_style=False)
            yaml_file.write("\n")

    # keep for all channels
    GLOBAL_NPs = []

    # loop first to extract the NPs
    RANKINGS = {}
    for channel in CHANNELS:

        RANKINGS[channel] = {}

        # folder
        if channel == "dplus":
            FOLDER = DPLUS_FOLDER
        elif channel == "dstar":
            FOLDER = DSTAR_FOLDER

        # fit results
        obs_fit_abs = "WCharm_lep_obs_OSSS_complete_xsec_alt_eta"
        f_result_abs = ROOT.TFile(os.path.join(FOLDER, obs_fit_abs, "Fits", f"{obs_fit_abs}.root"), "READ")
        f_result_abs_stat = ROOT.TFile(os.path.join(FOLDER, obs_fit_abs, "Fits", f"{obs_fit_abs}_statOnly.root"), "READ")
        fr_abs = f_result_abs.Get("nll_simPdf_newasimovData_with_constr")
        fr_abs_stat = f_result_abs_stat.Get("nll_simPdf_newasimovData_with_constr")

        # print
        for POI in POIs_abs:
            par = fr_abs.floatParsFinal().find(POI)
            par_stat = fr_abs_stat.floatParsFinal().find(POI)
            print(POI)
            print(par.getVal(), par.getErrorHi(), par.getErrorLo())
            print(par_stat.getVal(), par_stat.getErrorHi(), par_stat.getErrorLo())

        # read rankings
        # dict_keys(['Name', 'NPhat', 'NPerrHi', 'NPerrLo', 'POIup', 'POIdown', 'POIupPreFit', 'POIdownPreFit'])
        NPs = []
        for POI in POIs_abs:

            with open(os.path.join(FOLDER, obs_fit_abs, f"Ranking_{POI}.yaml"), 'r') as stream:
                ranking = yaml.safe_load(stream)
                RANKINGS[channel][POI] = ranking

                # calculate total error from ranking
                err_up_total = 0
                err_dn_total = 0
                for NP in ranking:
                    for key in ["POIup", "POIdown"]:
                        if NP[key] > 0:
                            err_up_total += NP[key]**2
                        else:
                            err_dn_total += NP[key]**2
                err_up_total = err_up_total**0.5
                err_dn_total = err_dn_total**0.5
                print("initial total error for: ", channel, POI, err_up_total, f"-{err_dn_total}")

                # keep NPs adding up to some threshold
                err_up = 0
                err_dn = 0
                threshold_reached = False
                for NP in ranking:
                    if not threshold_reached:
                        NPs += [NP["Name"]]
                        for key in ["POIup", "POIdown"]:
                            if NP[key] > 0:
                                err_up += NP[key]**2
                            else:
                                err_dn += NP[key]**2
                    if (err_up**0.5 / err_up_total) > PRUNING_THRESHOLD and (err_dn**0.5 / err_dn_total) > PRUNING_THRESHOLD:
                        threshold_reached = True

        # add to list of global NPs
        GLOBAL_NPs += NPs

    # loop again to make the tables
    for channel in CHANNELS:

        # folder
        if channel == "dplus":
            FOLDER = DPLUS_FOLDER
            br = 1.0 * 2.0
        elif channel == "dstar":
            FOLDER = DSTAR_FOLDER
            br = 0.677 * 2.0

        # fit results
        obs_fit_abs = "WCharm_lep_obs_OSSS_complete_xsec_alt_eta"
        f_result_abs = ROOT.TFile(os.path.join(FOLDER, obs_fit_abs, "Fits", f"{obs_fit_abs}.root"), "READ")
        f_result_abs_stat = ROOT.TFile(os.path.join(FOLDER, obs_fit_abs, "Fits", f"{obs_fit_abs}_statOnly.root"), "READ")
        fr_abs = f_result_abs.Get("nll_simPdf_newasimovData_with_constr")
        fr_abs_stat = f_result_abs_stat.Get("nll_simPdf_newasimovData_with_constr")

        # list of NPs (same for all channels)
        NPs = list(dict.fromkeys(GLOBAL_NPs))
        assert "MUON_SAGITTA_RESBIAS" in NPs, "MUON_SAGITTA_RESBIAS not included.. increase the threshold!"

        # get the remaining error of the pruned away NPs
        for POI in POIs_abs:

            # load the ranking
            ranking = RANKINGS[channel][POI]

            # keep NPs adding up to some threshold
            err_up_remaining = 0
            err_dn_remaining = 0
            count = 0
            for NP in ranking:
                if NP['Name'] not in NPs:
                    count += 1
                    # print(f"Adding up remaining error for {NP['Name']} for {POI}")
                    for key in ["POIup", "POIdown"]:
                        if NP[key] > 0:
                            err_up_remaining += NP[key]**2
                        else:
                            err_dn_remaining += NP[key]**2
            print("partial error for: ", channel, POI, err_up_remaining**0.5, f"-{err_dn_remaining**0.5}", count)
            if abs(err_up_remaining) > 0 or abs(err_dn_remaining) > 0:
                RANKINGS[channel][POI] += [{
                    'Name': f"Uncorr_{channel}_{POI}",
                    'POIup': err_up_remaining**0.5,
                    'POIdown': -err_dn_remaining**0.5,
                }]

        # NP impact table
        for charge in ["minus", "plus"]:

            # POIs
            if charge == "plus":
                POIs_channel = POIs_abs[0:5]
            elif charge == "minus":
                POIs_channel = POIs_abs[5:10]

            impact_table = {
                'independent_variables': [
                    {'header': {'name': 'LEP_ABS_ETA'},
                     'values': [{'high': bin_edges(POI, "eta")[0], 'low': bin_edges(POI, "eta")[1]} for POI in POIs_channel]}
                ],
                'dependent_variables': [
                    dependable_dict("DSIG/DLEP_ABS_ETA"),
                    dependable_dict("DSIG/DLEP_ABS_ETA (breakdown of systematics)"),
                ]
            }

            # fill in POIs
            for i, POI in enumerate(POIs_channel):
                par = fr_abs.floatParsFinal().find(POI)
                par_stat = fr_abs_stat.floatParsFinal().find(POI)
                vals_stat = atlas_rounding(par_stat.getVal() / br, par_stat.getErrorHi() / br, par_stat.getErrorLo() / br)
                vals = atlas_rounding(par.getVal() / br, (par.getErrorHi()**2 - par_stat.getErrorHi()**2)**0.5 / br,
                                      (par.getErrorLo()**2 - par_stat.getErrorLo()**2)**0.5 / br, vals_stat[-1])

                # create the two columns (2nd one gets filled out later)
                impact_table['dependent_variables'][0]['values'] += [{
                    'errors': [{'label': 'stat', 'symerror': vals_stat[1]},
                               {'label': 'sys', 'asymerror': {'minus': -vals[2], 'plus': vals[1]}}],
                    'value': vals[0]
                }]
                impact_table['dependent_variables'][1]['values'] += [{
                    'errors': [{'label': 'stat', 'symerror': vals_stat[1]}],
                    'value': vals[0]
                }]

            # sort NPs for each channel
            NPs_sorted = []
            for NP in NPs:
                sum_error = 0
                for POI in POIs_channel:
                    ranking = RANKINGS[channel][POI]
                    NP_impact = [x for x in ranking if x['Name'] == NP]
                    if len(NP_impact):
                        sum_error += max(abs(float(NP_impact[0]['POIup'])), abs(float(NP_impact[0]['POIdown'])))
                NPs_sorted += [(sum_error, NP)]
            NPs_sorted.sort(key=lambda x: x[0], reverse=True)

            # add the total correlated errors
            for POI in POIs_channel:
                NPs_sorted += [(0.0, f"Uncorr_{channel}_{POI}")]

            # fill rankings
            # total error validation
            err_up = {POI: 0 for POI in POIs_channel}
            err_dn = {POI: 0 for POI in POIs_channel}
            for _, NP in NPs_sorted:

                # loop over POIs
                for i, POI in enumerate(POIs_channel):

                    # load the ranking yaml
                    ranking = RANKINGS[channel][POI]

                    NP_impact = [x for x in ranking if x['Name'] == NP]
                    if len(NP_impact):
                        vals = atlas_rounding(par.getVal() / br, float(NP_impact[0]['POIup']) / br, float(NP_impact[0]['POIdown']) / br)
                    else:
                        # print(f"WARNING: NP {NP} not found for POI {POI}")
                        vals = 0, 0, 0
                    impact_table['dependent_variables'][1]['values'][i]['errors'] += [{'label': NP, 'asymerror': {'minus': vals[2], 'plus': vals[1]}}]

                    # validation
                    for i in [1, 2]:
                        if vals[i] > 0:
                            err_up[POI] += vals[i]**2
                        else:
                            err_dn[POI] += vals[i]**2

            for POI in POIs_channel:
                print(f"final total error for {channel} {POI}: {err_up[POI]**0.5} -{err_dn[POI]**0.5}")

            with open(f'hepdata/{channel}_{charge}.yaml', 'w') as yaml_file:
                yaml.dump(impact_table, yaml_file, default_flow_style=False)


if __name__ == "__main__":

    # run
    main()
