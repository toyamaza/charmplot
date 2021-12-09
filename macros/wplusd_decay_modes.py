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


def decode(i):
    if i == 1:
        return 111
    elif i == 2:
        return 211
    elif i == 3:
        return 321
    elif i == 4:
        return 310
    elif i == 5:
        return 130
    elif i == 6:
        return 11
    elif i == 7:
        return 12
    elif i == 8:
        return 13
    elif i == 9:
        return 14
    else:
        return 0


def main(options):

    # file
    f = ROOT.TFile(options.input, "READ")

    # histogram
    h = f.Get("MG_Wjets_emu_411MisMatched_OS-SS_0tag_Dplus_truth_Dmeson_full_decay_mode_dplus")
    total = h.GetSumOfWeights()
    decay_mode_dict = {}
    for i in range(0, h.GetNbinsX() + 2):
        y = h.GetBinContent(i)
        if y:
            decay_str = ""
            encoded = i
            while (encoded):
                num = encoded % 10
                decay_str += f"{decode(num)} "
                encoded = encoded // 10
            decay_mode_dict[decay_str] = y

    decay_mode_dict_sorted = dict(sorted(decay_mode_dict.items(), key=lambda item: item[1], reverse=True))
    for k, v in decay_mode_dict_sorted.items():
        print(f"{str(k):26}: {v:7.0f} ({100 * v / total:2.3f}%)")


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
