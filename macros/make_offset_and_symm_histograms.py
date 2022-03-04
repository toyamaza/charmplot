#!/usr/bin/env python
import ROOT
import re

def main(options, args):

    # histogram names
    dplus_histo_list = [
        "SS_el_minus_0tag_Dplus_Dmeson_m_fit",
        "SS_el_minus_1tag_Dplus_Dmeson_m_fit",
        "SS_el_plus_0tag_Dplus_Dmeson_m_fit",
        "SS_el_plus_1tag_Dplus_Dmeson_m_fit",
        "SS_mu_minus_0tag_Dplus_Dmeson_m_fit",
        "SS_mu_minus_1tag_Dplus_Dmeson_m_fit",
        "SS_mu_plus_0tag_Dplus_Dmeson_m_fit",
        "SS_mu_plus_1tag_Dplus_Dmeson_m_fit",
    ]
    dstar_histo_list = [
        "SS_el_minus_0tag_Dstar_Dmeson_mdiff_fit",
        "SS_el_minus_1tag_Dstar_Dmeson_mdiff_fit",
        "SS_el_plus_0tag_Dstar_Dmeson_mdiff_fit",
        "SS_el_plus_1tag_Dstar_Dmeson_mdiff_fit",
        "SS_mu_minus_0tag_Dstar_Dmeson_mdiff_fit",
        "SS_mu_minus_1tag_Dstar_Dmeson_mdiff_fit",
        "SS_mu_plus_0tag_Dstar_Dmeson_mdiff_fit",
        "SS_mu_plus_1tag_Dstar_Dmeson_mdiff_fit",
    ]

    pt_bin = ["pt_bin1", "pt_bin2", "pt_bin3", "pt_bin4", "pt_bin5"]

    data_list = []
    mc_list = []
    sig_list = []

    # Create histogram names for inclusive plots
    if (options.decay_mode == "Dplus"):
        data_list += ["Data_" + x for x in dplus_histo_list]
        mc_list += ["MC_TOT_" + x for x in dplus_histo_list]
        sig_list += ["Sherpa2211_WplusD_Matched_" + x for x in dplus_histo_list]
    else:
        data_list += ["Data_" + x for x in dstar_histo_list]
        mc_list += ["MC_TOT_" + x for x in dstar_histo_list]
        sig_list += ["Sherpa2211_WplusD_Matched_" + x for x in dstar_histo_list]

    # Create histogram names with pt bins included and add to data and mc list
    if (options.differential_bins):
        if (options.decay_mode == "Dplus"):
            for bin in pt_bin:
                data_list += ["Data_" + x.split("Dmeson")[0] + bin + x.split("Dplus")[-1] for x in dplus_histo_list]
                mc_list += ["MC_TOT_" + x.split("Dmeson")[0] + bin + x.split("Dplus")[-1] for x in dplus_histo_list]
                sig_list += ["Sherpa2211_WplusD_Matched_" + x.split("Dmeson")[0] + bin + x.split("Dplus")[-1] for x in dplus_histo_list]
        else:
            for bin in pt_bin:
                data_list += ["Data_" + x.split("Dmeson")[0] + bin + x.split("Dstar")[-1] for x in dstar_histo_list]
                mc_list += ["MC_TOT_" + x.split("Dmeson")[0] + bin + x.split("Dstar")[-1] for x in dstar_histo_list]
                sig_list += ["Sherpa2211_WplusD_Matched_" + x.split("Dmeson")[0] + bin + x.split("Dstar")[-1] for x in dstar_histo_list]

    # input file
    f = ROOT.TFile(options.input, "READ")

    # output file
    out_symm = ROOT.TFile("Mock_MC.root", "RECREATE")
    out_offset = ROOT.TFile("Offset.root", "RECREATE")
    out_sig_res = ROOT.TFile("Sig_Res.root", "RECREATE")

    # loop
    for data, mc, sig in zip(data_list, mc_list, sig_list):

        print(f"data = {data}")
        print(f"mc = {mc}")
        print(f"sig = {sig}")

        data_OS = data.replace("SS","OS")
        mc_OS = mc.replace("SS","OS")
        sig_OS = sig.replace("SS","OS")

        # Get the histogram
        f.cd()
        h_data = f.Get(data).Clone(data)
        h_data_OS = f.Get(data_OS).Clone(data_OS)
        h_mc = f.Get(mc).Clone(mc)
        h_mc_OS = f.Get(mc_OS).Clone(mc_OS)
        h_sig = f.Get(sig).Clone(sig)
        h_sig_OS = f.Get(sig_OS).Clone(sig_OS)

        # Create the Mock and Offset histograms
        h_symm_name = "Mock_MC_" + data.split("SS_")[-1]
        h_offset_name = "Offset_" + data.split("SS_")[-1]
        h_sig_bkg_diff_SS_name = "Sig_Res_" + data.split("SS_")[-1]
        h_sig_bkg_diff_OS_name = "Sig_Res_" + data.split("OS_")[-1]

        h_symm = h_data.Clone()
        h_symm.SetNameTitle(h_symm_name, h_symm_name)
        h_offset = h_data.Clone()
        h_offset.SetNameTitle(h_offset_name, h_offset_name)

        # Create Data - MC histogram
        h_diff = h_data.Clone(f"{h_data.GetName()}_diff")
        h_diff.Add(h_mc, -1)

        # Create Sig Histograms
        h_sig_bkg_diff_SS = h_diff.Clone()
        h_sig_bkg_diff_SS.SetNameTitle(h_sig_bkg_diff_SS_name, h_sig_bkg_diff_SS_name)
        h_sig_bkg_diff_SS.Add(h_sig,-1)
        h_sig_bkg_diff_OS = h_data_OS.Clone()
        h_sig_bkg_diff_OS.Add(h_mc_OS,-1)
        h_sig_bkg_diff_OS.Add(h_sig_OS,-1)
        h_sig_bkg_diff_OS.SetNameTitle(h_sig_bkg_diff_OS_name, h_sig_bkg_diff_OS_name)

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

        charge = re.findall("_minus_|_plus_", data)[0]
        flavor = re.findall("el_|mu_", data)[0][:-1]
        btag = re.findall("_0tag_|_1tag_", data)[0]
        meson = re.findall("Dplus_|Dstar_", data)[0]
        var = "Dmeson_" + data.split("Dmeson_")[-1].replace("_fit", "")
        ptbin = ""
        if "pt_bin" in data:
            ptbin = re.findall("_pt_bin[1-5]", data)[0]
        h_symm_OS_name = f"{flavor}{charge}SR{btag}{meson}OS_Mock_MC{ptbin}__{var}"
        h_symm_SS_name = f"{flavor}{charge}SR{btag}{meson}SS_Mock_MC{ptbin}__{var}"
        h_offset_OS_name = f"{flavor}{charge}SR{btag}{meson}OS_Offset{ptbin}__{var}"
        h_offset_SS_name = f"{flavor}{charge}SR{btag}{meson}SS_Offset{ptbin}__{var}"
        h_sig_res_OS_name = f"{flavor}{charge}SR{btag}{meson}OS_Sig_Res{ptbin}__{var}"
        h_sig_res_SS_name = f"{flavor}{charge}SR{btag}{meson}SS_Sig_Res{ptbin}__{var}"

        # Create OS and SS histograms
        h_symm_OS = h_symm.Clone()
        h_symm_OS.SetNameTitle(h_symm_OS_name, h_symm_OS_name)
        h_symm_SS = h_symm.Clone()
        h_symm_SS.SetNameTitle(h_symm_SS_name, h_symm_SS_name)
        h_offset_OS = h_offset.Clone()
        h_offset_OS.SetNameTitle(h_offset_OS_name, h_offset_OS_name)
        h_offset_SS = h_offset.Clone()
        h_offset_SS.SetNameTitle(h_offset_SS_name, h_offset_SS_name)
        h_sig_res_OS = h_sig_bkg_diff_OS.Clone()
        h_offset_SS.SetNameTitle(h_sig_res_OS_name, h_sig_res_OS_name)
        h_sig_res_SS = h_sig_bkg_diff_SS.Clone()
        h_offset_SS.SetNameTitle(h_sig_res_SS_name, h_sig_res_SS_name)

        # Save the Mock MC histogram
        out_symm.cd()
        h_symm_OS.Write()
        h_symm_SS.Write()

        # Save the offset histogram
        out_offset.cd()
        h_offset_OS.Write()
        h_offset_SS.Write()

        # Save the signal resolution histogram
        out_sig_res.cd()
        h_sig_res_OS.Write()
        h_sig_res_SS.Write()

    out_symm.Close()
    out_offset.Close()
    out_sig_res.Close()
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
    parser.add_option('--differential-bins',
                      action="store_true", dest="differential_bins",
                      default=False)

    # parse input arguments
    options, args = parser.parse_args()

    # run
    main(options, args)
