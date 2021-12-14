#!/usr/bin/env python
import ROOT
import re


def main(options, args):

    # histogram names
    dplus_histo_list = [
        "SS_el_minus_0tag_Dplus_Dmeson_m",
        "SS_el_minus_1tag_Dplus_Dmeson_m",
        "SS_el_plus_0tag_Dplus_Dmeson_m",
        "SS_el_plus_1tag_Dplus_Dmeson_m",
        "SS_mu_minus_0tag_Dplus_Dmeson_m",
        "SS_mu_minus_1tag_Dplus_Dmeson_m",
        "SS_mu_plus_0tag_Dplus_Dmeson_m",
        "SS_mu_plus_1tag_Dplus_Dmeson_m",
    ]
    dstar_histo_list = [
        "SS_el_minus_0tag_Dstar_Dmeson_mdiff",
        "SS_el_minus_1tag_Dstar_Dmeson_mdiff",
        "SS_el_plus_0tag_Dstar_Dmeson_mdiff",
        "SS_el_plus_1tag_Dstar_Dmeson_mdiff",
        "SS_mu_minus_0tag_Dstar_Dmeson_mdiff",
        "SS_mu_minus_1tag_Dstar_Dmeson_mdiff",
        "SS_mu_plus_0tag_Dstar_Dmeson_mdiff",
        "SS_mu_plus_1tag_Dstar_Dmeson_mdiff",
    ]

    pt_bin = ["pt_bin1", "pt_bin2", "pt_bin3", "pt_bin4", "pt_bin5"]

    data_list = []
    mc_list = []

    # Create histogram names for inclusive plots
    if (options.decay_mode == "Dplus"):
        data_list += ["Data_" + x for x in dplus_histo_list]
        mc_list += ["MC_TOT_" + x for x in dplus_histo_list]
    else:
        data_list += ["Data_" + x for x in dstar_histo_list]
        mc_list += ["MC_TOT_" + x for x in dstar_histo_list]

    # Create histogram names with pt bins included and add to data and mc list
    if (options.decay_mode == "Dplus"):
        for bin in pt_bin:
            data_list += ["Data_" + x.split("Dmeson")[0] + bin + x.split("Dplus")[-1] for x in dplus_histo_list]
            mc_list += ["MC_TOT_" + x.split("Dmeson")[0] + bin + x.split("Dplus")[-1] for x in dplus_histo_list]
    else:
        for bin in pt_bin:
            data_list += ["Data_" + x.split("Dmeson")[0] + bin + x.split("Dstar")[-1] for x in dstar_histo_list]
            mc_list += ["MC_TOT_" + x.split("Dmeson")[0] + bin + x.split("Dstar")[-1] for x in dstar_histo_list]

    # input file
    f = ROOT.TFile(options.input, "READ")

    # output file
    out_symm = ROOT.TFile("Mock_MC.root", "RECREATE")
    out_offset = ROOT.TFile("Offset.root", "RECREATE")

    # loop
    for data, mc in zip(data_list, mc_list):

        # Get the histogram
        f.cd()
        h_data = f.Get(data).Clone(data)
        h_mc = f.Get(mc).Clone(mc)

        # Create the Mock and Offset histograms
        h_symm_name = "Mock_MC_" + data.split("SS_")[-1]
        h_offset_name = "Offset_" + data.split("SS_")[-1]

        h_symm = h_data.Clone()
        h_symm.SetNameTitle(h_symm_name, h_symm_name)
        h_offset = h_data.Clone()
        h_offset.SetNameTitle(h_offset_name, h_offset_name)

        # Create Data - MC histogram
        h_diff = h_data.Clone(f"{h_data.GetName()}_diff")
        h_diff.Add(h_mc, -1)

        # Fill histograms
        for i in range(1, h_data.GetNbinsX() + 1):

            # Set Data-MC and data uncertainty tolerance (set to 20 x data uncertainty, 20 sigma)
            data_mc_diff = h_diff.GetBinContent(i)
            data_unc_tol = 20 * h_data.GetBinError(i)

            if (data_mc_diff > data_unc_tol):
                h_symm.SetBinContent(i, data_mc_diff)
                h_offset.SetBinContent(i, 1e-4)
            elif (data_mc_diff < data_unc_tol and data_mc_diff > 0.0):
                h_symm.SetBinContent(i, data_unc_tol)
                h_offset.SetBinContent(i, (data_unc_tol - abs(data_mc_diff)))
            else:
                h_symm.SetBinContent(i, data_unc_tol)
                h_offset.SetBinContent(i, (data_unc_tol + abs(data_mc_diff)))

            # Set symmetric histogram error
            h_symm.SetBinError(i, h_diff.GetBinError(i))

        if re.search("pt_bin", data):
            h_symm_OS_name = re.split("[0, 1]tag", data.split("SS_")[-1])[0] + "SR" + re.split("(minus|plus)", data.split("_pt_bin")[0])[-1] + "_OS_Mock_MC" + re.split("(Dplus|Dstar)", data.split("_Dmeson")[0])[-1] + "_" + re.split("_pt_bin[1-5]", data)[-1]
            h_symm_SS_name = re.split("[0, 1]tag", data.split("SS_")[-1])[0] + "SR" + re.split("(minus|plus)", data.split("_pt_bin")[0])[-1] + "_SS_Mock_MC" + re.split("(Dplus|Dstar)", data.split("_Dmeson")[0])[-1] + "_" + re.split("_pt_bin[1-5]", data)[-1]
            h_offset_OS_name = re.split("[0, 1]tag", data.split("SS_")[-1])[0] + "SR" + re.split("(minus|plus)", data.split("_pt_bin")[0])[-1] + "_OS_Offset" + re.split("(Dplus|Dstar)", data.split("_Dmeson")[0])[-1] + "_" + re.split("_pt_bin[1-5]", data)[-1]
            h_offset_SS_name = re.split("[0, 1]tag", data.split("SS_")[-1])[0] + "SR" + re.split("(minus|plus)", data.split("_pt_bin")[0])[-1] + "_SS_Offset" + re.split("(Dplus|Dstar)", data.split("_Dmeson")[0])[-1] + "_" + re.split("_pt_bin[1-5]", data)[-1]
        else:
            h_symm_OS_name = re.split("[0, 1]tag", data.split("SS_")[-1])[0] + "SR" + re.split("(minus|plus)", data.split("_Dmeson")[0])[-1] + "_OS_Mock_MC__Dmeson_" + data.split("_Dmeson_")[-1]
            h_symm_SS_name = re.split("[0, 1]tag", data.split("SS_")[-1])[0] + "SR" + re.split("(minus|plus)", data.split("_Dmeson")[0])[-1] + "_SS_Mock_MC__Dmeson_" + data.split("_Dmeson_")[-1]
            h_offset_OS_name = re.split("[0, 1]tag", data.split("SS_")[-1])[0] + "SR" + re.split("(minus|plus)", data.split("_Dmeson")[0])[-1] + "_OS_Offset__Dmeson_" + data.split("_Dmeson_")[-1]
            h_offset_SS_name = re.split("[0, 1]tag", data.split("SS_")[-1])[0] + "SR" + re.split("(minus|plus)", data.split("_Dmeson")[0])[-1] + "_SS_Offset__Dmeson_" + data.split("_Dmeson_")[-1]

        # Create OS and SS histograms
        h_symm_OS = h_symm.Clone()
        h_symm_OS.SetNameTitle(h_symm_OS_name, h_symm_OS_name)
        h_symm_SS = h_symm.Clone()
        h_symm_SS.SetNameTitle(h_symm_SS_name, h_symm_SS_name)
        h_offset_OS = h_offset.Clone()
        h_offset_OS.SetNameTitle(h_offset_OS_name, h_offset_OS_name)
        h_offset_SS = h_offset.Clone()
        h_offset_SS.SetNameTitle(h_offset_SS_name, h_offset_SS_name)

        # Save the Mock MC histogram
        out_symm.cd()
        h_symm_OS.Write()
        h_symm_SS.Write()

        # Save the offset histogram
        out_offset.cd()
        h_offset_OS.Write()
        h_offset_SS.Write()

    out_symm.Close()
    out_offset.Close()
    f.Close()


if __name__ == "__main__":

    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-i', '--input',
                      action="store", dest="input",
                      help="Path to the prefit histograms.root file")
    parser.add_option('-d', '--decay',
                      action="store", dest="decay_mode",
                      help="Decay Mode, ie. Dplus, Dstar")

    # parse input arguments
    options, args = parser.parse_args()

    # run
    main(options, args)
