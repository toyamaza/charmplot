#!/usr/bin/env python
import ROOT

# file names
backgrounds = {
    "Wjets_emu_Rest": {
        "inputs": ["MG_Wjets_cjets_emu", "MG_Wjets_bjets_emu", "MG_Wjets_light_emu"],
        "regions": ["HardMisMatched", "Other"],
    },
    "Wjets_emu_MisMatched": {
        "inputs": ["MG_Wjets_cjets_emu", "MG_Wjets_bjets_emu", "MG_Wjets_light_emu"],
        "regions": ["MisMatched", "MatchedNoFid"],
    },
    "Other": {
        "inputs": ["Diboson", "MG_Zjets_light_emu", "MG_Zjets_cjets_emu", "MG_Zjets_bjets_emu", "MG_Zjets_tau", "MG_Wjets_tau"],
        "regions": [""],
    },
    "MultiJet": {
        "inputs": ["data"],
        "regions": [""],
        "matrix_method": True,
    },
}

# multijet
matrix_method = {}

# scale factors
channels = ["el_minus_SR_0tag_Dplus", "el_plus_SR_0tag_Dplus", "mu_minus_SR_0tag_Dplus", "mu_plus_SR_0tag_Dplus",
            "el_minus_Anti_SR_0tag_Dplus", "el_plus_Anti_SR_0tag_Dplus", "mu_minus_Anti_SR_0tag_Dplus", "mu_plus_Anti_SR_0tag_Dplus"]

# pt bins
ptbins = ['pt_bin1', 'pt_bin2', 'pt_bin3', 'pt_bin4', 'pt_bin5', '']

# var
var = "Dmeson_m"

# output
f_out = ROOT.TFile("inclusive_bkg.root", "RECREATE")
f_out.Close()

# loop
for bkg, props in backgrounds.items():
    for ptbin in ptbins:
        for charge in ["OS", "SS"]:

            # output file
            f_out = ROOT.TFile("inclusive_bkg.root", "UPDATE")

            # input files
            files = [ROOT.TFile(f"{f}.root", "READ") for f in props["inputs"]]
            print(files)

            # histogram
            h = None
            for f in files:
                for c in channels:
                    for r in props["regions"]:
                        h_name = f"{c}_{charge}{('_' + r if r else '')}{('_' + ptbin if ptbin else '')}__{var}"
                        print(f"reading histogram {h_name} from file {f}...")
                        if "matrix_method" not in props:
                            h_temp = f.Get(h_name)
                        else:
                            h_temp = None
                            h_temp_Tight = f.Get(f"Tight_{h_name}")
                            h_temp_AntiTight = f.Get(f"AntiTight_{h_name}")
                            if h_temp_Tight:
                                h_temp = h_temp_Tight.Clone(h_temp_Tight.GetName().replace("Tight_", ""))
                                if h_temp_AntiTight:
                                    h_temp.Add(h_temp_AntiTight)
                        if h_temp:
                            print("got histogram")
                            if not h:
                                h = h_temp.Clone(f"{bkg}_{charge}{('_' + ptbin if ptbin else '')}__{var}")
                            else:
                                h.Add(h_temp)
                        else:
                            print(f"WARNING:: histogram {h_name} not found in file {f}!!")

            # write
            f_out.cd()
            h.Write()
            f_out.Close()


f_out.Close()
