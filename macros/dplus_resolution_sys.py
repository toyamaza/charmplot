#!/usr/bin/env python
import os
import ROOT
import shutil
dirname = os.path.join(os.path.dirname(__file__), "../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasLabels.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasUtils.C"))
ROOT.SetAtlasStyle()

rand = ROOT.TRandom3()

parameters = ROOT.TFile("/global/cfs/cdirs/atlas/wcharm/charmpp_data/systematics/track_eff_systematics.v3.root", "READ")


def createCanvasPads(name):
    c = ROOT.TCanvas(name, name, 1200, 1200)
    # Upper histogram plot is pad1
    pad1 = ROOT.TPad(f"pad1_{name}", f"pad1_{name}", 0, 0.30, 1, 1.0)
    pad1.SetTopMargin(0.1)
    pad1.SetBottomMargin(0.05 / 0.70)
    pad1.SetFillStyle(4000)
    pad1.Draw()
    # Lower ratio plot is pad2
    c.cd()
    pad2 = ROOT.TPad(f"pad2_{name}", f"pad2_{name}", 0, 0.00, 1, 0.40)
    pad2.SetTopMargin(0.05 / 0.40)
    pad2.SetBottomMargin(0.10 / 0.40)
    pad2.SetFillStyle(4000)
    pad2.Draw()

    return c, pad1, pad2


def configure_axis(h_up, h_dn):
    GLOBAL_SF = 1.0
    h_up.GetYaxis().SetTitleSize(h_up.GetYaxis().GetTitleSize() * GLOBAL_SF)
    h_up.GetYaxis().SetTitleOffset(h_up.GetYaxis().GetTitleOffset() * (1 / (GLOBAL_SF * 1.1)))
    h_up.GetYaxis().SetLabelSize(h_up.GetYaxis().GetLabelSize() * GLOBAL_SF)
    h_up.GetXaxis().SetLabelSize(0)

    SF = 0.70 / 0.40
    h_dn.GetYaxis().SetTitleSize(h_dn.GetYaxis().GetTitleSize() * SF * GLOBAL_SF)
    h_dn.GetYaxis().SetTitleOffset(h_dn.GetYaxis().GetTitleOffset() * (1 / (GLOBAL_SF * 1.06 * SF)))
    h_dn.GetXaxis().SetTitleSize(h_dn.GetXaxis().GetTitleSize() * SF * GLOBAL_SF)
    h_dn.GetXaxis().SetTitleOffset(h_dn.GetXaxis().GetTitleOffset() * (1 / (GLOBAL_SF * 0.6 * SF)))
    h_dn.GetXaxis().SetLabelOffset(h_dn.GetXaxis().GetLabelOffset() * 4.0)
    h_dn.GetYaxis().SetLabelSize(h_dn.GetYaxis().GetLabelSize() * SF * GLOBAL_SF)
    h_dn.GetXaxis().SetLabelSize(h_dn.GetXaxis().GetLabelSize() * SF * GLOBAL_SF * 1.2)

    # tick marks
    h_up.GetYaxis().SetNdivisions(506)
    h_dn.GetYaxis().SetNdivisions(306)


def atlas_label(labels=[]):
    l1 = ROOT.TLatex()
    l1.SetTextFont(73)
    l1.SetTextSize(42)
    l1.DrawLatexNDC(0.18, 0.82, "ATLAS")
    l2 = ROOT.TLatex()
    l2.SetTextFont(43)
    l2.SetTextSize(42)
    l2.DrawLatexNDC(0.33, 0.82, "Internal")
    for i, lab in enumerate(labels):
        l2.DrawLatexNDC(0.18, 0.82 - (i + 1) * 0.05, lab)


def make_legend(N):
    leg = ROOT.TLegend(0.58, 0.848 - N / 2. * (0.095), 0.92, 0.848)
    leg.SetBorderSize(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(42)
    leg.SetTextFont(43)
    return leg


def main(options):

    if "Dplus" in options.decay:
        shutil.copy("Sh_WplusD.root", "Sh_WplusD.root.BAK")
        shutil.copy("Sh_WplusD.root", "Sh_WplusD_2.root")

        f = ROOT.TFile("Sh_WplusD_2.root", "UPDATE")

        var = "m"
    else:
        shutil.copy("MG_NLO_WplusD.root", "MG_NLO_WplusD.root.BAK")
        shutil.copy("MG_NLO_WplusD.root", "MG_NLO_WplusD_2.root")

        f = ROOT.TFile("MG_NLO_WplusD_2.root", "UPDATE")

        var = "mdiff"

    for i in range(1, 6):
        for sign in ["OS", "SS"]:
            for lep in ["el", "mu"]:
                for charge in ["minus", "plus"]:
                    for pt_bin in [f"pt_bin{j}" for j in range(1, 6)]:

                        f.cd()
                        name = f"{lep}_{charge}_SR_0tag_{options.decay}_{sign}_Matched_truth_pt_bin{i}_{pt_bin}__Dmeson_{var}"
                        h = f.Get(name)
                        print(f"histogram name = {name}")
                        if not h:
                            continue

                        h_sys = f.Get(f"DMESON_RESO_ADHOC_4_-_{name}")
                        h_sys.SetLineColor(ROOT.kRed)

                        h_tot = f.Get(f"DMESON_RESO_TOT_-_{name}")
                        h_tot.SetLineColor(ROOT.kBlue)

                        if "Dplus" in options.decay:
                            if charge == "minus":
                                resolution = parameters.Get(f"pars_DPlus_nominal_{pt_bin}").GetBinContent(6)
                            else:
                                resolution = parameters.Get(f"pars_DMinus_nominal_{pt_bin}").GetBinContent(6)

                            smear = parameters.Get(f"Dplus_TRACK_EFF_TOT_{pt_bin}_sigma_diff").getVal() / 1000.
                        else:
                            if charge == "minus":
                                resolution = parameters.Get(f"pars_DstarPlus_nominal_{pt_bin}").GetBinContent(6)
                            else:
                                resolution = parameters.Get(f"pars_DstarMinus_nominal_{pt_bin}").GetBinContent(6)

                            smear = parameters.Get(f"Dstar_TRACK_EFF_TOT_{pt_bin}_sigma_diff").getVal() / 1000.

                        print(f"Smearing {options.decay} {sign} {lep} {charge} {pt_bin}: {resolution} {smear}")

                        # ad-hoc uncertainty
                        for sigma in [1, 2, 4]:
                            h_adhoc = h.Clone(f"RESO_ADHOC_{sigma}_-_{name}")
                            for k in range(1, h_adhoc.GetNbinsX() + 1):
                                h_adhoc.SetBinContent(k, 0)

                            kernels = []
                            for k in range(1, h_adhoc.GetNbinsX() + 1):
                                h_temp = h_adhoc.Clone(f"{h.GetName()}_{k}")
                                m = h_adhoc.GetXaxis().GetBinCenter(k)
                                for _ in range(10000):
                                    diff = rand.Gaus(0.0, resolution / sigma)
                                    m_smear = m + diff
                                    h_temp.Fill(m_smear)
                                h_temp.Scale(h.GetBinContent(k) / h_temp.Integral())
                                kernels += [h_temp]

                            for h_temp in kernels:
                                h_adhoc.Add(h_temp)

                            # write new histogram
                            h_adhoc.Write(f"RESO_ADHOC_{sigma}_-_{name}")
                            print(f"Wrote {h_adhoc.GetName()}")

                        # track resolution uncertainty
                        h_tot_new = h.Clone(f"RESO_TOT_-_{name}")
                        for k in range(1, h_tot_new.GetNbinsX() + 1):
                            h_tot_new.SetBinContent(k, 0)

                        kernels = []
                        for k in range(1, h_tot_new.GetNbinsX() + 1):
                            h_temp = h_tot_new.Clone(f"{h.GetName()}_{k}")
                            m = h_tot_new.GetXaxis().GetBinCenter(k)
                            for _ in range(10000):
                                diff = rand.Gaus(0.0, resolution / sigma)
                                m_smear = m + diff
                                h_temp.Fill(m_smear)
                            h_temp.Scale(h.GetBinContent(k) / h_temp.Integral())
                            kernels += [h_temp]

                        for h_temp in kernels:
                            h_tot_new.Add(h_temp)

                        # write new histogram
                        h_tot_new.Write(f"RESO_TOT_-_{name}")
                        print(f"Wrote {h_tot_new.GetName()}")

                        configure_axis(h, h_sys)

                        # legend
                        leg = make_legend(3)
                        leg.AddEntry(h, "Nominal", "f")
                        leg.AddEntry(h_tot, "Total material sys", "f")
                        leg.AddEntry(h_sys, "ad-hoc reso (1/4 #sigma)", "f")

                        h_sys_ratio = h_sys.Clone()
                        h_sys_ratio.Divide(h)

                        h_tot_ratio = h_tot.Clone()
                        h_tot_ratio.Divide(h)

                        c, pad1, pad2 = createCanvasPads(f"{name}")
                        pad1.cd()
                        h.Draw("hist")
                        h_sys.Draw("hist same")
                        h_tot.Draw("hist same")
                        h.SetMaximum(1.5 * h.GetMaximum())

                        atlas_label([f"{sign} 0tag {lep} {charge}", f"truth_pt_bin{i}"])
                        leg.Draw()

                        pad2.cd()
                        pad2.SetGridy()
                        h_sys_ratio.Draw("hist")
                        h_tot_ratio.Draw("hist same")
                        h_sys_ratio.GetYaxis().SetRangeUser(0.81, 1.19)

                        if "Dplus" in options.decay:
                            c.Print(f"dplus_reso/{name}_truth_pt_bin{i}.pdf")
                        else:
                            c.Print(f"dstar_reso/{name}_truth_pt_bin{i}.pdf")


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-i', '--input',
                      action="store", dest="input")
    parser.add_option('-d', '--decay',
                      action="store", dest="decay")

    # parse input arguments
    options, _ = parser.parse_args()

    # make output dir
    if "Dplus" in options.decay:
        if not os.path.isdir("dplus_reso"):
            os.makedirs("dplus_reso")
    if "Dstar" in options.decay:
        if not os.path.isdir("dstar_reso"):
            os.makedirs("dstar_reso")

    # run
    main(options)
