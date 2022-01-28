#!/usr/bin/env python
import os
import ROOT

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasLabels.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasUtils.C"))
ROOT.SetAtlasStyle()


def decay_modes(pdgId, index):
    if (pdgId == 411):
        if index == 1:
            return "Other", 0.0
        elif index == 2:
            return "Kππ", 0.0938
        elif index == 3:
            return "KKπ", 0.00968
        elif index == 4:
            return "πππ", 0.00327
        elif index == 5:
            return "Kπππ0", 0.0625
        elif index == 6:
            return "ππππ0", 0.0116
        elif index == 7:
            return "Kπμν", 0.0365
        elif index == 8:
            return "Kπeν", 0.0402
        elif index == 9:
            return "K0πππ", 0.062
    elif (pdgId == 421):
        if index == 1:
            return "Other", 0.0
        elif index == 2:
            return "Kπ", 0.03950
        elif index == 3:
            return "Kππ0", 0.144
        elif index == 4:
            return "Kππ0π0", 0.0886
        elif index == 5:
            return "K0ππ", 0.056
        elif index == 6:
            return "Kπππ", 0.0823
        elif index == 7:
            return "Kππππ0", 0.043
        elif index == 8:
            return "ππππ", 0.00756
        elif index == 9:
            return "πππππ0", 0.0042
    elif (pdgId == 431):
        if index == 1:
            return "Other", 0.0
        elif index == 2:
            return "KKπ", 0.0539
        elif index == 3:
            return "Kππ", 0.0065
        elif index == 4:
            return "πππ", 0.0108
        elif index == 5:
            return "πππɣ", 0.0394 * (0.425 * 0.3941 + 0.295) + 0.168 * 0.0422
        elif index == 6:
                            # η(π+π-π0)π+     η(π+π-ɣ)π+π0               η'(π+π-ɣ/π+π-η(ɣɣ))π+π0
            return "ππππ0", 0.0168 * 0.2292 + (0.089 + 0.095) * 0.0422 + (0.058 + 0.056) * (0.295 + 0.425 * 0.3941)
        elif index == 7:
            return "πππX", 0
        elif index == 8:
            return "KKμν", 0.0116871
        elif index == 9:
            return "KKeν", 0.0116871
    return "", 0.0


def main(options):

    # file
    f = ROOT.TFile(options.input, "READ")

    # histograms
    h_dplus_MG = f.Get("MG_Wjets_Dmeson_D_decay_dplus")
    h_dplus_Sh = f.Get("Sherpa2211_Wjets_Dmeson_D_decay_dplus")
    h_dzero_MG = f.Get("MG_Wjets_Dmeson_D_decay_dzero")
    h_dzero_Sh = f.Get("Sherpa2211_Wjets_Dmeson_D_decay_dzero")
    h_dsubs_MG = f.Get("MG_Wjets_Dmeson_D_decay_dsubs")
    h_dsubs_Sh = f.Get("Sherpa2211_Wjets_Dmeson_D_decay_dsubs")

    print("~~~ D+ ~~~")
    print(f"{'':10s}{'Pythia':10s}{'Sherpa':10s}{'PDG':10s}")
    for i in range(1, h_dplus_MG.GetNbinsX() + 1):
        print(f"{decay_modes(411, i)[0]:10s}{h_dplus_MG.GetBinContent(i) / h_dplus_MG.GetSumOfWeights():1.4f}{'':4s}{h_dplus_Sh.GetBinContent(i) / h_dplus_Sh.GetSumOfWeights():1.4f}{'':4s}{decay_modes(411, i)[1]:1.5f}")

    print("~~~ D0 ~~~")
    print(f"{'':10s}{'Pythia':10s}{'Sherpa':10s}{'PDG':10s}")
    for i in range(1, h_dzero_MG.GetNbinsX() + 1):
        print(f"{decay_modes(421, i)[0]:10s}{h_dzero_MG.GetBinContent(i) / h_dzero_MG.GetSumOfWeights():1.4f}{'':4s}{h_dzero_Sh.GetBinContent(i) / h_dzero_Sh.GetSumOfWeights():1.4f}{'':4s}{decay_modes(421, i)[1]:1.5f}")

    print("~~~ Ds ~~~")
    print(f"{'':10s}{'Pythia':10s}{'Sherpa':10s}{'PDG':10s}")
    for i in range(1, h_dsubs_MG.GetNbinsX() + 1):
        print(f"{decay_modes(431, i)[0]:10s}{h_dsubs_MG.GetBinContent(i) / h_dsubs_MG.GetSumOfWeights():1.4f}{'':4s}{h_dsubs_Sh.GetBinContent(i) / h_dsubs_Sh.GetSumOfWeights():1.4f}{'':4s}{decay_modes(431, i)[1]:1.5f}")


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-i', '--input',
                      action="store", dest="input",
                      default="", help="Input file.")

    # parse input arguments
    options, _ = parser.parse_args()

    # do the plotting
    main(options)
