#!/usr/bin/env python
import ROOT


def main(options, args):
    # file names
    backgrounds = {
        "MG_Wjets_Rest": {
            "inputs": ["MG_Wjets_cjets_emu", "MG_Wjets_bjets_emu", "MG_Wjets_light_emu", "MG_Wjets_tau"],
            "regions": ["HardMisMatched", "Other"],
        },
        "MG_Wjets_MisMatched": {
            "inputs": ["MG_Wjets_cjets_emu", "MG_Wjets_bjets_emu", "MG_Wjets_light_emu", "MG_Wjets_tau"],
            "regions": ["MisMatched", "MatchedNoFid"],
        },
        "Sh_Wjets_Rest": {
            "inputs": ["Sh_2211_Wjets_light_emu", "Sh_2211_Wjets_cjets_emu", "Sh_2211_Wjets_bjets_emu", "Sh_2211_Wjets_tau"],
            "regions": ["HardMisMatched", "Other"],
        },
        "Sh_Wjets_MisMatched": {
            "inputs": ["Sh_2211_Wjets_light_emu", "Sh_2211_Wjets_cjets_emu", "Sh_2211_Wjets_bjets_emu", "Sh_2211_Wjets_tau"],
            "regions": ["MisMatched", "MatchedNoFid"],
        },
    }

    # scale factors
    channels = [f"el_minus_SR_0tag_{options.decay_mode}", f"el_plus_SR_0tag_{options.decay_mode}", f"mu_minus_SR_0tag_{options.decay_mode}", f"mu_plus_SR_0tag_{options.decay_mode}"]
    channels_loose = [f"el_minus_Anti_SR_0tag_{options.decay_mode}", f"el_plus_Anti_SR_0tag_{options.decay_mode}", f"mu_minus_Anti_SR_0tag_{options.decay_mode}", f"mu_plus_Anti_SR_0tag_{options.decay_mode}"]

    # pt bins
    ptbins = ['pt_bin1', 'pt_bin2', 'pt_bin3', 'pt_bin4', 'pt_bin5', '']

    # var
    if (options.decay_mode == "Dplus"):
        var = "Dmeson_m"
    else:
        var = "Dmeson_mdiff"

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

                # correct scaling
                sf = None
                if len(props["regions"]) > 1:
                    norm = {x: 0 for x in props["regions"]}
                    norm_loose = {x: 0 for x in props["regions"]}
                    for r in props["regions"]:
                        for f in files:
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
                h_loose = None
                for r in props["regions"]:
                    for f in files:
                        for c in channels + channels_loose:
                            h_name = f"{c}_{charge}{('_' + r if r else '')}{('_' + ptbin if ptbin else '')}__{var}"
                            print(f"reading histogram {h_name} from file {f}...")
                            if "matrix_method" not in props:
                                h_temp = f.Get(h_name)
                                if h_temp and sf and c in channels_loose:
                                    print(f"scaling histogram by {sf[r]}")
                                    h_temp.Scale(sf[r])
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
                                if "Anti" in c:
                                    if not h_loose:
                                        h_loose = h_temp.Clone(f"{bkg}_{charge}_loose{('_' + ptbin if ptbin else '')}__{var}")
                                    else:
                                        h_loose.Add(h_temp)
                            else:
                                print(f"WARNING:: histogram {h_name} not found in file {f}!!")

                # write
                f_out.cd()
                h.Write()
                h_loose.Write()
                f_out.Close()

    f_out.Close()


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-d', '--decay-mode',
                      action="store", dest="decay_mode",
                      default="", help="Decay Mode: Dplus, Dstar, etc.")

    # parse input arguments
    options, args = parser.parse_args()

    # run
    main(options, args)
