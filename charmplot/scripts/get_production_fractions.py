#!/usr/bin/env python
import os
import ROOT

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.SetAtlasStyle()


if __name__ == "__main__":

    samples = ["sherpa", "madgraph", "powheg"]

    species = [
        "D",
        "D0",
        "DS",
        "LambdaC",
        "Cascade",
        "Cascade0",
        "OmegaC",
    ]

    files = {}
    files["sherpa"] = ROOT.TFile("Sh_Wjets.root")
    files["madgraph"] = ROOT.TFile("MG_Wjets.root")
    files["powheg"] = ROOT.TFile("Ph_Wjets.root")

    production = {s: {} for s in samples}
    production_from_b = {s: {} for s in samples}

    for s in samples:
        f = files[s]
        f.cd("truth")
        for c in species:
            production[s][c] = 0
            production_from_b[s][c] = 0
            h = f.Get(f"truth/truth__N_{c}")
            hb = f.Get(f"truth/truth__N_{c}_from_B")
            for i in range(1, h.GetNbinsX() + 1):
                production[s][c] += h.GetBinContent(i) * (i - 1)
                production_from_b[s][c] += hb.GetBinContent(i) * (i - 1)

    total = {s: 0 for s in samples}
    for s in samples:
        for c in species:
            total[s] += production[s][c]
            total[s] += production_from_b[s][c]

    print("==== non-B production =====")

    txt = "\t".join([s for s in samples])
    print(f"\t{txt}")
    for c in species:
        txt = "\t".join(["%.4e" % (production[s][c]) for s in samples])
        print(f"{c}\t{txt}")

    print("==== production from B =====")

    txt = "\t".join([s for s in samples])
    print(f"\t{txt}")
    for c in species:
        txt = "\t".join(["%.4e" % (production_from_b[s][c]) for s in samples])
        print(f"{c}\t{txt}")

    print("==== relative non-B production =====")

    txt = "\t".join([s for s in samples])
    print(f"\t{txt}")
    for c in species:
        txt = "\t".join(["%.4e" % (100 * production[s][c] / total[s]) for s in samples])
        print(f"{c}\t{txt}")

    print("==== relative production from B =====")

    txt = "\t".join([s for s in samples])
    print(f"\t{txt}")
    for c in species:
        txt = "\t".join(["%.4e" % (100 * production_from_b[s][c] / total[s]) for s in samples])
        print(f"{c}\t{txt}")
