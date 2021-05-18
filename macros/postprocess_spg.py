#!/usr/bin/env python
import re
import ROOT

# file names
names = [
    "SPG_Baryonbar",
    "SPG_Baryon",
    "SPG_Dsminus",
    "SPG_Ds",
    "SPG_D0bar",
    "SPG_D0",
]

# scale factors
scale_factors = {
    "Ds_OS": [0.3599, 0.3046, 0.3313, 0.3989, 0.4803],
    "Ds_SS": [0.0840, 0.0389, 0.0313, 0.0329, 0.0503],
    "D0_OS": [0.2725, 0.1940, 0.1782, 0.2045, 0.2620],
    "D0_SS": [0.1804, 0.1588, 0.1667, 0.1776, 0.1867],
    "Baryon_OS": [0.3287, 0.2000, 0.2398, 0.2615, 0.3204],
    "Baryon_SS": [0.1193, 0.0384, 0.0226, 0.0235, 0.0354],
}

# vars
accepted_vars = ["Dmeson_m", "Dmeson_mdiff", "Dmeson_pt"]

# loop
for name in names:
    flavor = name.replace("SPG_", "").replace("bar", "").replace("minus", "")
    f_in = ROOT.TFile(f"{name}.root", "READ")
    f_out = ROOT.TFile(f"{name}_postProc.root", "RECREATE")

    # do OS-SS
    subtract = False
    if flavor in ["Ds", "Baryon"]:
        subtract = True

    # construct inclusive histograms from scaled ones
    inclusive_hists = {}

    # loop through histograms
    for key in f_in.GetListOfKeys():

        # TObject name
        obj_name = key.GetName()

        # var
        var = obj_name.split("__")[-1]
        if var not in accepted_vars:
            continue

        # pt bin
        pt_bin = re.findall("pt_bin([0-9])", obj_name)

        # if not inclusive
        if pt_bin:

            # should be length 1
            pt_bin = int(pt_bin[0])

            # charge (OS or SS)
            charge = re.findall("[Dplustar]+_([OS]+)", obj_name)[0]

            # print out
            print(obj_name, pt_bin, charge, flavor, scale_factors[f"{flavor}_{charge}"][pt_bin - 1])

            # scale and save
            if not subtract:
                h_scaled = f_in.Get(obj_name)
                h_scaled.Scale(scale_factors[f"{flavor}_{charge}"][pt_bin - 1])
                f_out.cd()
                h_scaled.Write(h_scaled.GetName())
            else:
                if not (f_in.Get(obj_name.replace("SS", "OS")) and f_in.Get(obj_name.replace("OS", "SS"))):
                    continue
                h_OS = f_in.Get(obj_name.replace("SS", "OS")).Clone(f"{obj_name}_OS")
                h_SS = f_in.Get(obj_name.replace("OS", "SS")).Clone(f"{obj_name}_SS")
                h_scaled = h_OS.Clone(f"{obj_name}_OS-SS")
                h_scaled.Add(h_SS, -1.0)
                print(f"created {h_scaled} from {h_OS} - {h_SS}")
                h_scaled.Scale(scale_factors[f"{flavor}_{charge}"][pt_bin - 1])
                f_out.cd()
                h_scaled.Write(obj_name)

            # inclusive name
            inc_name = obj_name.replace(f"_pt_bin{pt_bin}", "")
            if inc_name not in inclusive_hists:
                h_inc = h_scaled.Clone(inc_name)
                inclusive_hists[inc_name] = h_inc
                print(f"{inc_name} created...")
            else:
                inclusive_hists[inc_name].Add(h_scaled)
                print(f"added {h_scaled} to {inc_name}...")

    # save inclusive histograms
    for name, h in inclusive_hists.items():
        f_out.cd()
        h.Write(name)

    f_in.Close()
    f_out.Close()
