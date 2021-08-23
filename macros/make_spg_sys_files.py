#!/usr/bin/env python
import ROOT

# file names
backgrounds = {
    "Wjets_emu_Charm": {
        "inputs": [
            ["MG_Wjets_emu"],
        ],
        "regions": [
            ["411MisMatched", "413MisMatched", "421MisMatched", "431MisMatched", "BaryonMisMatched"],
        ]
    },
}

# nominal spg
spg = ["SPG_Dminus_postProc", "SPG_Dplus_postProc", "SPG_Ds_postProc", "SPG_D0_Dstar_postProc",
       "SPG_Baryon_postProc", "SPG_Dsminus_postProc", "SPG_D0_Dstar_bar_postProc", "SPG_Baryonbar_postProc"]

# spg sys
systematics = {
    "Dplus": {"SPG_Dplus_postProc": 1.05, "SPG_Dminus_postProc": 1.05},
    "Dzero": {"SPG_D0_Dstar_postProc": 1.05, "SPG_D0_Dstar_bar_postProc": 1.05},
    "Dsubs": {"SPG_Ds_postProc": 1.05, "SPG_Dsminus_postProc": 1.05},
}

# scale factors
channels = ["el_minus_SR_0tag_Dplus", "el_plus_SR_0tag_Dplus", "mu_minus_SR_0tag_Dplus", "mu_plus_SR_0tag_Dplus"]

# pt bins
ptbins = ['pt_bin1', 'pt_bin2', 'pt_bin3', 'pt_bin4', 'pt_bin5', '']

# var
var = "Dmeson_m"

# output
f_out = ROOT.TFile("spg_sys.root", "RECREATE")
f_out.Close()

# sys output
for syst in systematics:
    f_out = ROOT.TFile(f"spg_sys_{syst}.root", "RECREATE")
    f_out.Close()

# loop
for bkg, props in backgrounds.items():
    for ptbin in ptbins:
        for charge in ["OS", "SS"]:

            # spg file
            h_spg = None
            h_spg_sys = {key: None for key in systematics.keys()}
            spg_files = []
            for file in spg:
                f = ROOT.TFile(f"{file}.root", "READ")
                spg_files += [f]
                h_name = f"inclusive_Dplus_{charge}{('_' + ptbin if ptbin else '')}__{var}"
                h_temp = f.Get(h_name)
                if h_temp:
                    if not h_spg:
                        h_spg = h_temp.Clone(h_name + "_tmp")
                    else:
                        h_spg.Add(h_temp)
                for syst in systematics:
                    if h_temp:
                        h_temp_sys = h_temp.Clone(f"{h_temp.GetName()}_tmp1")
                        for k, v in systematics[syst].items():
                            if file == k:
                                h_temp_sys.Scale(v)
                                print(f"Scaling sys histogram {h_temp_sys} by {v}")
                                break
                        if not h_spg_sys[syst]:
                            h_spg_sys[syst] = h_temp_sys.Clone(h_name + f"_tmp_{syst}")
                        else:
                            h_spg_sys[syst].Add(h_temp_sys)

            for c in channels:

                # histogram
                h = None
                all_files = []
                for i in range(len(props["inputs"])):
                    for r in props["regions"][i]:
                        for file in props["inputs"][i]:
                            f = ROOT.TFile(f"{file}.root", "READ")
                            all_files += [f]
                            h_name = f"{c}_{charge}{('_' + r if r else '')}{('_' + ptbin if ptbin else '')}__{var}"
                            # print(f"reading histogram {h_name} from file {f}...")
                            h_temp = f.Get(h_name)
                            if h_temp:
                                h_temp = h_temp.Clone(f"{h_temp.GetName()}_clone")
                            if h_temp:
                                # print(f"got histogram {h_temp}")
                                if not h:
                                    h = h_temp.Clone(f"{c}_{charge}_{bkg}{('_' + ptbin if ptbin else '')}__{var}")
                                else:
                                    h.Add(h_temp)
                            else:
                                # print(f"WARNING:: histogram {h_name} not found in file {f}!!")
                                pass

                # write
                f_out = ROOT.TFile("spg_sys.root", "UPDATE")
                f_out.cd()
                h_spg_out = h_spg.Clone(h_spg.GetName() + "_tmp")
                SF = h.GetSumOfWeights() / h_spg_out.GetSumOfWeights()
                h_spg_out.Scale(SF)
                h_spg_out.Write(h.GetName())
                print(f"SAVED {h.GetName()} integral: {h_spg_out.GetSumOfWeights()}")
                f_out.Close()
                for syst in systematics:
                    f_out = ROOT.TFile(f"spg_sys_{syst}.root", "UPDATE")
                    f_out.cd()
                    h_spg_out_sys = h_spg_sys[syst].Clone(h_spg.GetName() + "_tmp_" + syst)
                    h_spg_out_sys.Scale(SF)
                    h_spg_out_sys.Write(h.GetName())
                    print(f"SAVED {h.GetName()} integral: {h_spg_out_sys.GetSumOfWeights()} for sys {syst}")
                    f_out.Close()

                _ = [f.Close() for f in all_files]

            _ = [f.Close() for f in spg_files]


f_out.Close()
