#!/usr/bin/env python
import os
import ROOT
import yaml

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasLabels.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasUtils.C"))
ROOT.SetAtlasStyle()

DPLUS_FOLDER = "/global/cscratch1/sd/mmuskinj/TRExFitter/fit_2022_12_14_Dplus_new"
DSTAR_FOLDER = "/global/cscratch1/sd/mmuskinj/TRExFitter/fit_2022_12_14_Dstar_new"
# DPLUS_FOLDER = "/global/cfs/cdirs/atlas/wcharm/TRExFitter/Output/Dplus_2022_08_05_v2"
# DSTAR_FOLDER = "/global/cfs/cdirs/atlas/wcharm/TRExFitter/Output/Dstar_2022_08_11_v2"
# DPLUS_FOLDER = "/global/cfs/cdirs/atlas/wcharm/TRExFitter/Output/Dplus_2022_07_26_fullRanking_v2"
# DSTAR_FOLDER = "/global/cfs/cdirs/atlas/wcharm/TRExFitter/Output/Dstar_2022_08_08_fullRanking_v2"

CHANNELS = ["dplus", "dstar"]

ETA_BINS = [f"W^{{-}} bin {i}" for i in range(1, 6)] + [f"W^{{+}} bin {i}" for i in range(1, 6)]

# make output folder
if not os.path.isdir("covariance"):
    os.makedirs("covariance")


def main():

    # read inputs
    for channel in CHANNELS:

        # pt eta
        for var in ["pt", "eta"]:

            # folder
            if channel == "dplus":
                FOLDER = DPLUS_FOLDER
                meson = "Dplus"
            elif channel == "dstar":
                FOLDER = DSTAR_FOLDER
                meson = "Dstar"
            obs_fit_abs = f"WCharm_{meson}_lep_obs_OSSS_complete_xsec_alt_{var}"

            # for cov in ["Covariance_statOnly", "Combined_Covariance_postfit"]:
            for cov in ["Covariance_statOnly", "Combined_Covariance_prefit", "Combined_Covariance_postfit"]:
                with open(os.path.join(FOLDER, obs_fit_abs, "Covariances", f"{cov}.yaml"), 'r') as stream:
                    cov_stat = yaml.safe_load(stream)
                    print(cov_stat)
                    h_cov = ROOT.TH2D(f"h_{cov}_{channel}_{var}", f"h_{cov}_{channel}_{var}", 10, 0, 10, 10, 0, 10)
                    h_cov.SetMarkerSize(0.8)
                    h_cov.GetZaxis().SetTitle("Covariance [pb^{2}]")
                    h_cov.GetZaxis().SetTitleSize(h_cov.GetZaxis().GetTitleSize() * 0.8)
                    h_cov.GetZaxis().SetTitleOffset(h_cov.GetZaxis().GetTitleOffset() * 1.5)
                    h_cov.GetZaxis().SetLabelSize(h_cov.GetZaxis().GetLabelSize() * 0.8)
                    if var == "eta":
                        h_cov.GetXaxis().SetTitle("|#eta(lep)| bin")
                        h_cov.GetYaxis().SetTitle("|#eta(lep)| bin")
                    elif var == "pt":
                        h_cov.GetXaxis().SetTitle("p_{T}(D) bin")
                        h_cov.GetYaxis().SetTitle("p_{T}(D) bin")
                    h_cov.GetXaxis().SetTitleSize(h_cov.GetXaxis().GetTitleSize() * 0.8)
                    h_cov.GetXaxis().SetTitleOffset(h_cov.GetXaxis().GetTitleOffset() * 1.2)
                    h_cov.GetXaxis().SetLabelOffset(h_cov.GetXaxis().GetLabelOffset() * 1.2)
                    h_cov.GetYaxis().SetTitleSize(h_cov.GetYaxis().GetTitleSize() * 0.8)
                    h_cov.GetYaxis().SetTitleOffset(h_cov.GetYaxis().GetTitleOffset() * 1.4)
                    for i in range(10):
                        h_cov.GetXaxis().SetBinLabel(i + 1, ETA_BINS[i])
                        h_cov.GetYaxis().SetBinLabel(i + 1, ETA_BINS[i])
                        for j in range(10):
                            if channel == "dplus":
                                h_cov.SetBinContent(i + 1, j + 1, cov_stat[1]["covariance_rows"][i][j] / 4.)
                            elif channel == "dstar":
                                h_cov.SetBinContent(i + 1, j + 1, cov_stat[1]["covariance_rows"][i][j] / (2 * 0.677)**2)

                    ROOT.gStyle.SetPaintTextFormat("1.2e")
                    ROOT.gStyle.SetPalette(ROOT.kCool)
                    ROOT.gStyle.SetNumberContours(256)
                    canv = ROOT.TCanvas(f"{cov}_{channel}_{var}", f"{cov}_{channel}_{var}", 1000, 800)
                    canv.SetRightMargin(0.18)
                    canv.SetLeftMargin(0.18)
                    canv.SetBottomMargin(0.20)
                    canv.SetTopMargin(0.10)
                    h_cov.Draw("COLZ text")
                    if cov != "Covariance_statOnly":
                        # h_cov.SetMinimum(0)
                        # h_cov.SetMaximum(0.65)
                        # h_cov.GetZaxis().SetRangeUser(0.0, 0.7)
                        pass
                    else:
                        # h_cov.GetZaxis().SetRangeUser(-0.001, 0.035)
                        pass

                    # ATLAS label
                    ROOT.ATLASLabel(0.20, 0.95, "", 1, 0.035)
                    if cov == "Covariance_statOnly":
                        ROOT.myText(0.50, 0.95, 1, "Stat. Only Covariance", 0.035)
                    elif cov == "Combined_Covariance_prefit":
                        ROOT.myText(0.50, 0.95, 1, "Combined Pre-fit Covariance", 0.035)
                    elif cov == "Combined_Covariance_postfit":
                        ROOT.myText(0.50, 0.95, 1, "Combined Post-fit Covariance", 0.035)
                    ROOT.myText(0.20, 0.91, 1, "#sqrt{s} = 13 TeV, 140 fb^{-1}", 0.035)
                    if channel == "dplus":
                        ROOT.myText(0.50, 0.91, 1, "W(#rightarrowl#nu)+D(#rightarrowK#pi#pi)", 0.035)
                    elif channel == "dstar":
                        ROOT.myText(0.50, 0.91, 1, "W(#rightarrowl#nu)+D*(#rightarrow(K#pi)#pi)", 0.035)

                    print(f" ================== {cov}_{channel}_{var} ==================")
                    cov_string = """dependent_variables:
  - header: {name: COV, units: 'pb$^{2}$'}
    qualifiers:
      - name: SQRT(s)
        units: GeV
        value: 13000
      - name: LUMINOSITY
        units: fb$^{-1}$
        value: 139
    values:\n"""
                    for i in range(1, 11):
                        for j in range(1, 11):
                            cov_string += f"      - value: {h_cov.GetBinContent(i, j):.2e}\n"
                    cov_string += """independent_variables:
  - header:
      name: bin1
    values:\n"""

                    for i in range(1, 11):
                        for j in range(1, 11):
                            cov_string += f"      - value: {h_cov.GetXaxis().GetBinLabel(i)}\n"
                    cov_string += """  - header:
      name: bin2
    values:\n"""
                    for i in range(1, 11):
                        for j in range(1, 11):
                            cov_string += f"      - value: {h_cov.GetXaxis().GetBinLabel(j)}\n"

                    # save output file
                    with open(f"covariance/{cov}_{channel}_{var}.yaml", 'w') as stream:
                        stream.write(cov_string)
                    canv.Print(f"covariance/{cov}_{channel}_{var}.pdf")


if __name__ == "__main__":

    # run
    main()
