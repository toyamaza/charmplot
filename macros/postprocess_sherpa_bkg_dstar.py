#!/usr/bin/env python
import ROOT

# file names
backgrounds = {
    "Wjets_emu_Bkg": {
        "inputs": [
            ["Sh_Wjets_light_emu", "Sh_Wjets_highpt_emu", "Sh_Wjets_cjets_emu", "Sh_Wjets_bjets_emu"],
            ["MG_Wjets_emu"],
        ],
        "regions": [
            ["HardMisMatched"],
            ["MisMatched", "MatchedNoFid", "Other"],
        ]
    },
    "Wjets_emu_Rest": {
        "inputs": [
            ["Sh_Wjets_light_emu", "Sh_Wjets_highpt_emu", "Sh_Wjets_cjets_emu", "Sh_Wjets_bjets_emu"],
            ["MG_Wjets_emu"],
        ],
        "regions": [
            ["HardMisMatched"],
            ["Other"],
        ]
    },
}

# multijet
matrix_method = {}

# scale factors
channels = ["el_minus_SR_0tag_Dstar", "el_plus_SR_0tag_Dstar", "mu_minus_SR_0tag_Dstar", "mu_plus_SR_0tag_Dstar"]
channels_loose = ["el_minus_Anti_SR_0tag_Dstar", "el_plus_Anti_SR_0tag_Dstar", "mu_minus_Anti_SR_0tag_Dstar", "mu_plus_Anti_SR_0tag_Dstar"]

# pt bins
ptbins = ['pt_bin1', 'pt_bin2', 'pt_bin3', 'pt_bin4', 'pt_bin5', '']

# var
var = "Dmeson_mdiff"

# output
f_out = ROOT.TFile("inclusive_sherpa_bkg.root", "RECREATE")
f_out.Close()

# loop
for bkg, props in backgrounds.items():
    for ptbin in ptbins:
        for charge in ["OS", "SS"]:

            # correct scaling
            sf = None
            if len(props["regions"]) > 1:
                regions = [item for sublist in props["regions"] for item in sublist]
                norm = {x: 0 for x in regions}
                norm_loose = {x: 0 for x in regions}
                for i in range(len(props["inputs"])):
                    files = []
                    for r in props["regions"][i]:
                        for file in props["inputs"][i]:
                            f = ROOT.TFile(f"{file}.root", "READ")
                            files += [f]
                            for c in channels:
                                h_name = f"{c}_{charge}{('_' + r if r else '')}{('_' + ptbin if ptbin else '')}__{var}"
                                h_temp = f.Get(h_name)
                                if h_temp:
                                    norm[r] += h_temp.GetSumOfWeights()
                            for c in channels_loose:
                                h_name = f"{c}_{charge}{('_' + r if r else '')}{('_' + ptbin if ptbin else '')}__{var}"
                                h_temp = f.Get(h_name)
                                if h_temp:
                                    norm_loose[r] += h_temp.GetSumOfWeights()
                    _ = [f.Close() for f in files]

                norm_sum = sum([norm[x] for x in norm])
                norm_sum_loose = sum([norm_loose[x] for x in norm_loose])
                sf = {x: (norm[x] / norm_sum) / (norm_loose[x] / norm_sum_loose) for x in norm}
                print("--- scaling between loose and SR ---")
                print(f"{charge} {ptbin} {bkg}")
                print(norm)
                print(norm_loose)
                print(sf)

            # histogram
            h = None
            all_files = []
            for i in range(len(props["inputs"])):
                for r in props["regions"][i]:
                    for file in props["inputs"][i]:
                        f = ROOT.TFile(f"{file}.root", "READ")
                        all_files += [f]
                        for c in channels + channels_loose:
                            h_name = f"{c}_{charge}{('_' + r if r else '')}{('_' + ptbin if ptbin else '')}__{var}"
                            print(f"reading histogram {h_name} from file {f}...")
                            h_temp = f.Get(h_name)
                            if h_temp:
                                h_temp = h_temp.Clone(f"{h_temp.GetName()}_clone")
                            if h_temp and sf and c in channels_loose:
                                print(f"scaling histogram by {sf[r]}")
                                h_temp.Scale(sf[r])
                            if h_temp:
                                print(f"got histogram {h_temp}")
                                if not h:
                                    h = h_temp.Clone(f"{bkg}_{charge}{('_' + ptbin if ptbin else '')}__{var}")
                                else:
                                    h.Add(h_temp)
                            else:
                                print(f"WARNING:: histogram {h_name} not found in file {f}!!")

            # write
            f_out = ROOT.TFile("inclusive_sherpa_bkg.root", "UPDATE")
            f_out.cd()
            h.Write()
            f_out.Close()
            _ = [f.Close() for f in all_files]

            # save per-channel histograms
            print("~~~~ saving per-channel histograms ~~~~")
            for c in channels:
                all_files = []
                for i in range(len(props["inputs"])):
                    for r in props["regions"][i]:
                        # output file
                        h = None
                        for file in props["inputs"][i]:
                            f = ROOT.TFile(f"{file}.root", "READ")
                            all_files += [f]
                            h_name = f"{c}_{charge}{('_' + r if r else '')}{('_' + ptbin if ptbin else '')}__{var}"
                            print(f"reading histogram {h_name} from file {f}...")
                            h_temp = f.Get(h_name)
                            if h_temp:
                                h_temp = h_temp.Clone(f"{h_temp.GetName()}_clone")
                                print(f"got histogram {h_temp} {h_temp.GetSumOfWeights()}")
                                if not h:
                                    h = h_temp.Clone(f"{c}_{charge}_{r}{('_' + ptbin if ptbin else '')}__{var}")
                                else:
                                    h.Add(h_temp)
                            else:
                                print(f"WARNING:: histogram {h_name} not found in file {f}!!")

                        # write
                        f_out = ROOT.TFile("inclusive_sherpa_bkg.root", "UPDATE")
                        f_out.cd()
                        print(f"Written histogram {h}")
                        # Check if histogram exists, skip otherwise
                        if h:
                            h.Write()
                        f_out.Close()
                _ = [f.Close() for f in all_files]


f_out.Close()
