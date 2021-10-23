#!/usr/bin/env python
from charmplot.common import utils
from charmplot.control import globalConfig
from charmplot.control import sample
import copy
import os
import ROOT

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasLabels.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasUtils.C"))
ROOT.SetAtlasStyle()

# samples
sample_c = sample.Sample("Wjets_cjets_ftag", None, **{'add': [], 'subtract': [], 'legendLabel': "true c-jets (rw)", 'lineColor': "ROOT.kGreen+1"})
sample_b = sample.Sample("Wjets_bjets_ftag", None, **{'add': [], 'subtract': [], 'legendLabel': "true b-jets", 'lineColor': "ROOT.kRed"})
sample_l = sample.Sample("Wjets_light_ftag", None, **{'add': [], 'subtract': [], 'legendLabel': "true l-jets", 'lineColor': "ROOT.kBlue"})
sample_c_dplus = sample.Sample("Wjets_cjets_dplus_ftag", None, **{'add': [], 'subtract': [], 'legendLabel': "true c-jets (D+)", 'lineColor': "ROOT.kMagenta+1"})
sample_c_dzero = sample.Sample("Wjets_cjets_dzero_ftag", None, **{'add': [], 'subtract': [], 'legendLabel': "true c-jets (D0)", 'lineColor': "ROOT.kMagenta+2"})
sample_c_dsubs = sample.Sample("Wjets_cjets_dsubs_ftag", None, **{'add': [], 'subtract': [], 'legendLabel': "true c-jets (Ds)", 'lineColor': "ROOT.kMagenta+3"})
sample_c_baryon = sample.Sample("Wjets_cjets_baryon_ftag", None, **{'add': [], 'subtract': [],
                                'legendLabel': "true c-jets (#Lambda)", 'lineColor': "ROOT.kBlue"})
sample_c_norw = sample.Sample("PROD_FRAC_NOWEIGHT_Wjets_cjets_ftag", None, **
                              {'add': [], 'subtract': [], 'legendLabel': "true c-jets", 'lineColor': "ROOT.kRed+2"})
all_samples = [
    sample_c,
    sample_b,
    # sample_l,
    sample_c_dplus,
    sample_c_dzero,
    sample_c_dsubs,
    sample_c_baryon,
    sample_c_norw,
]


def main(options, conf):

    # fitted variable
    var = conf.get_var(options.var)
    print(var.name)

    # channels
    channel_notag = conf.get_channel(f"{options.channel}_notag_ftag")
    channel_tag = conf.get_channel(f"{options.channel}_tag_ftag")
    print(channel_notag.name, channel_tag.name)

    # label
    channel_tag.label = ["DL1r, btag@70%", "jet p_{T} > 20 GeV, |#eta| < 2.5"]

    # file
    files = {}
    for inpt in options.input.split(","):
        f = ROOT.TFile(os.path.join(inpt, "histograms.root"), "READ")
        files[inpt] = f

    # get histos
    proxy = None
    histos = {inpt: {"notag": {}, "tag": {}} for inpt in options.input.split(",")}
    for inpt in options.input.split(","):
        for s in [x.name for x in all_samples]:
            hadron = ""
            if "dplus" in s or "dzero" in s or "dsubs" in s or "baryon" in s:
                hadron = "_hadron"
            print(f"{s}_{channel_notag.name}{hadron}_{var.name}")
            h_notag = files[inpt].Get(f"{s}_{channel_notag.name}{hadron}_{var.name}")
            h_tag = files[inpt].Get(f"{s}_{channel_tag.name}{hadron}_{var.name}")
            histos[inpt]["notag"][s] = h_notag
            histos[inpt]["tag"][s] = h_tag
            if proxy is None:
                proxy = h_notag.Clone("proxy")

    # calculate efficiency
    mc_map = {inpt: {} for inpt in options.input.split(",")}
    samples = [sample_b, sample_c_dplus, sample_c_dzero, sample_c_dsubs, sample_c_baryon, sample_c, sample_c_norw]
    for inpt in options.input.split(","):
        for s in samples:
            h_eff = histos[inpt]["tag"][s.name].Clone(f'{histos[inpt]["tag"][s.name].GetName()}_{os.path.basename(inpt).replace("yaml", "")}_eff')
            h_eff.Divide(histos[inpt]["notag"][s.name])
            for i in range(1, h_eff.GetNbinsX() + 1):
                y = h_eff.GetBinContent(i)
                e = h_eff.GetBinError(i)
                h_eff.SetBinContent(i, y / (y + 1))
                h_eff.SetBinError(i, e / (y + 1)**2)

            # mc map
            h_eff.SetLineColor(s.lineColor)
            h_eff.SetMarkerSize(0)
            mc_map[inpt][s] = h_eff

    samples = [sample_b, sample_c_dplus, sample_c_dzero, sample_c_dsubs, sample_c_baryon, sample_c]
    for inpt in options.input.split(","):
        # plot efficiency
        canv = utils.make_canvas_mc_ratio(proxy, var, channel_tag, "Ratio", x=800, y=800, events="b-tag efficiency",
                                          y_split=0, suffix=os.path.basename(inpt).replace("yaml", ""))
        canv.pad1.cd()

        # make legend
        canv.make_legend(mc_map[inpt], samples, print_yields=False)

        # set maximum after creating legend
        canv.set_maximum([mc_map[inpt][s] for s in samples], var, mc_map[inpt][samples[0]])

        proxy.GetYaxis().SetRangeUser(0, 1.0)
        for s in mc_map[inpt]:
            mc_map[inpt][s].Draw("hist e same")

        canv.print(f"{inpt}/efficiency.pdf")

    # combined samples (total)
    NAMES = ["individual", "total"]
    SAMPLES = [
        [sample_c_dplus, sample_c_dzero, sample_c_dsubs, sample_c_baryon],
        [sample_c, sample_c_norw],
    ]
    for name, samples in zip(NAMES, SAMPLES):
        mc_map_combined = {}
        samples_comb = []
        for s in samples:
            for inpt in options.input.split(","):
                samples_comb += [copy.deepcopy(s)]
                s_new = samples_comb[-1]
                mc_map_combined[s_new] = mc_map[inpt][s]
                # s_new.legendLabel = s_new.legendLabel.replace("true c-jets", "c")
                if "sherpa" in inpt:
                    s_new.legendLabel += " Sh2.2.1"
                elif "madgraph" in inpt:
                    s_new.legendLabel += " MGPy8"
                if "sherpa" in inpt:
                    mc_map_combined[s_new].SetLineStyle(2)
                    mc_map_combined[s_new].SetLineWidth(4)

        # plot efficiency combined
        canv = utils.make_canvas_mc_ratio(proxy, var, channel_tag, "Ratio", x=800, y=800, events="b-tag efficiency", y_split=0.3, suffix=f"{name}_total")
        canv.pad1.cd()

        # make legend
        canv.make_legend(mc_map_combined, samples_comb, print_yields=False, leg_offset=-0.10)

        # set maximum after creating legend
        canv.set_maximum([mc_map_combined[s] for s in samples_comb], var, mc_map_combined[samples_comb[0]])

        proxy.GetYaxis().SetRangeUser(0, 1.0)
        for s in mc_map_combined:
            mc_map_combined[s].Draw("hist e same")

        # ratio plot
        h_tmp = []
        canv.pad2.cd()
        if len(samples_comb) == 4:
            canv.set_ratio_range(0.61, 1.19, override=True)
        denum = mc_map_combined[samples_comb[0]]
        for i in range(0, len(mc_map_combined)):
            h_eff_ratio = mc_map_combined[samples_comb[i]].Clone(f"{mc_map_combined[samples_comb[i]].GetName()}_ratio")
            h_eff_ratio.Divide(denum)
            h_tmp += [h_eff_ratio]
            h_eff_ratio.Draw("hist same")
        # for i in range(0, len(mc_map_combined), 2):
        #     h_eff_ratio = mc_map_combined[samples_comb[i]].Clone(f"{mc_map_combined[samples_comb[i]].GetName()}_ratio")
        #     h_eff_ratio.Divide(mc_map_combined[samples_comb[i + 1]])
        #     h_tmp += [h_eff_ratio]
        #     h_eff_ratio.Draw("hist e same")

        canv.print(f"{inpt}/efficiency_comb_{name}.pdf")

    # close
    for _, f in files.items():
        f.Close()


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-a', '--analysis-config',
                      action="store", dest="analysis_config",
                      help="analysis config file")
    parser.add_option('-i', '--input',
                      action="store", dest="input",
                      default="", help="Input folder path.")
    parser.add_option('-c', '--channel',
                      action="store", dest="channel",
                      help="channel to use",
                      default="lep_1jet")
    parser.add_option('-v', '--var',
                      action="store", dest="var",
                      help="fitted variable",
                      default="jet_pt")

    # parse input arguments
    options, args = parser.parse_args()

    # analysis configs
    config = options.analysis_config

    # config object
    conf = globalConfig.GlobalConfig(config, "")

    # do the plotting
    main(options, conf)
