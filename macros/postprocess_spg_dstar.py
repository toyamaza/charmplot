#!/usr/bin/env python
import os
import re
import ROOT
import sys

# file names
names = [
    "SPG_Dplus",
    "SPG_Dminus",
    "SPG_Baryonbar",
    "SPG_Baryon",
    "SPG_Dsminus",
    "SPG_Ds",
    "SPG_D0bar",
    "SPG_D0",
    "SPG_Dstar_plus",
    "SPG_Dstar_minus",
]

# name map
name_map = {
    "D": "411MisMatched",
    "Ds": "431MisMatched",
    "D0": "421MisMatched",
    "Baryon": "BaryonMisMatched",
    "Dstar": "413MisMatched"
}

# keep only the following truth categories
truth_cat_to_keep = {
    "Dstar": "413MisMatched",
}

# vars
# accepted_vars = ["Dmeson_m", "Dmeson_mdiff", "Dmeson_pt"]
accepted_vars = ["Dmeson_mdiff"]


def main(options, args):

    if options.scale_factors and os.path.isfile(options.scale_factors):
        f_scale_factors = ROOT.TFile(options.scale_factors, "READ")
        print(f"Got the scale factors file {f_scale_factors}")
    else:
        print("No scale factors provided. SPG will not be scaled.")

    for name in names:
        flavor = name.replace("SPG_", "").replace("_bar", "").replace("bar", "").replace("minus", "").replace("plus", "").replace("Dstar_","Dstar")
        f_in = ROOT.TFile(f"{name}.root", "READ")
        f_out = ROOT.TFile(f"{name}_postProc.root", "RECREATE")

        # do OS-SS
        OS_only = False
        if flavor in ["D", "Ds", "Baryon"]:
            OS_only = True

        # construct inclusive histograms from scaled ones
        inclusive_hists = {}

        # loop through histograms
        for key in f_in.GetListOfKeys():

            # TObject name
            obj_name = key.GetName()

            # check truth_cat_to_keep
            if flavor in truth_cat_to_keep:
                if truth_cat_to_keep[flavor] not in obj_name:
                    continue

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
                decay_mode = re.findall("([Dplustar]+)_[OS]+", obj_name)[0]
                charge = re.findall("[Dplustar]+_([OS]+)", obj_name)[0]

                # var for scaling
                if decay_mode == "Dplus":
                    scale_var = "Dmeson_m"
                elif decay_mode == "Dstar":
                    scale_var = "Dmeson_mdiff"

                # print out
                print(obj_name, pt_bin, charge, flavor)

                # scale factors
                sf = 1.0
                if options.scale_factors:
                    print(f"Getting scale-factors for {flavor}_{charge} in pt_bin{pt_bin}")
                    name_spg = f"SPG_{name_map[flavor]}_{charge}_{decay_mode}_{name_map[flavor]}_pt_bin{pt_bin}_{scale_var}"
                    name_mg = f"Wjets_emu_{name_map[flavor]}_{charge}_{decay_mode}_{name_map[flavor]}_pt_bin{pt_bin}_{scale_var}"
                    h_mg = f_scale_factors.Get(name_mg)
                    h_spg = f_scale_factors.Get(name_spg)
                    if not (h_spg and h_mg):
                        print(f"Unable to retreive histograms {name_spg} and {name_mg} from file {f_scale_factors}! Exiting!")
                        sys.exit(1)
                    sf = h_mg.GetSumOfWeights() / h_spg.GetSumOfWeights()

                # scale and save
                if not OS_only:
                    h_scaled = f_in.Get(obj_name)
                    if options.scale_factors:
                        h_scaled.Scale(sf)
                    f_out.cd()
                else:
                    if not (f_in.Get(obj_name.replace("SS", "OS")) and f_in.Get(obj_name.replace("OS", "SS"))):
                        continue
                    h_OS = f_in.Get(obj_name.replace("SS", "OS")).Clone(f"{obj_name}_OS")
                    # h_SS = f_in.Get(obj_name.replace("OS", "SS")).Clone(f"{obj_name}_SS")
                    h_scaled = h_OS.Clone(f"{obj_name}_OS_only")
                    # h_scaled.Add(h_SS, -1.0)
                    print(f"created {h_scaled} from {h_OS}")
                    if options.scale_factors:
                        print(f"sf: {sf}")
                        h_scaled.Scale(sf)
                    f_out.cd()

                # save to file
                h_scaled.Write(obj_name)
                if flavor in truth_cat_to_keep:
                    h_scaled.Write(obj_name.replace(f"_{truth_cat_to_keep[flavor]}", ""))

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
            if flavor in truth_cat_to_keep:
                h.Write(name.replace(f"_{truth_cat_to_keep[flavor]}", ""))

        f_in.Close()
        f_out.Close()

    if options.scale_factors:
        f_scale_factors.Close()


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-s', '--scale-factors',
                      action="store", dest="scale_factors",
                      default="", help="Path to the histograms.root file from the initial spg comparison."
                      "If empty, no scaling will be performed")

    # parse input arguments
    options, args = parser.parse_args()

    # run
    main(options, args)
