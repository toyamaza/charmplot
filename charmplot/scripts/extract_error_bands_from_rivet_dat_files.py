#!/usr/bin/env python
from charmplot.common import utils
from charmplot.control import channel
from charmplot.control import sample
from charmplot.control import variable
import array
import logging
import os
import re
import ROOT
import sys

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.SetAtlasStyle()

# logging
root = logging.getLogger()
root.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
root.addHandler(handler)

# proxy samples
samples_dict = {
    "PDF260000": sample.Sample('PDF260000', None, **{'add': [], 'subtract': [], 'legendLabel': 'NNPDF30_nlo_as_0118', 'lineColor': 'ROOT.kBlack'}),
    "PDF261000": sample.Sample('PDF261000', None, **{'add': [], 'subtract': [], 'legendLabel': 'NNPDF30_nnlo_as_0118', 'lineColor': 'ROOT.kGreen+4'}),
    "PDF303600": sample.Sample('PDF303600', None, **{'add': [], 'subtract': [], 'legendLabel': 'NNPDF31_nnlo_as_0118', 'lineColor': 'ROOT.kGreen+2'}),
    "PDF25200": sample.Sample('PDF25200', None, **{'add': [], 'subtract': [], 'legendLabel': 'MMHT2014nlo68clas118', 'lineColor': 'ROOT.kGreen+0'}),
    "PDF65060": sample.Sample('PDF65060', None, **{'add': [], 'subtract': [], 'legendLabel': 'ATLAS-epWZ16-EIG', 'lineColor': 'ROOT.kYellow+4'}),
    "PDF14400": sample.Sample('PDF14400', None, **{'add': [], 'subtract': [], 'legendLabel': 'CT18NLO', 'lineColor': 'ROOT.kYellow+2'}),
    "PDF14600": sample.Sample('PDF14600', None, **{'add': [], 'subtract': [], 'legendLabel': 'CT18ANLO', 'lineColor': 'ROOT.kYellow+0'}),
    "QCD": sample.Sample('QCD', None, **{'add': [], 'subtract': [], 'legendLabel': 'QCDScale', 'lineColor': 'ROOT.kBlue'}),
    "LOMG": sample.Sample('LOMG', None, **{'add': [], 'subtract': [], 'legendLabel': 'LO MG', 'lineColor': 'ROOT.kRed'}),
    # "MGFxFx": sample.Sample('MGFxFx', None, **{'add': [], 'subtract': [], 'legendLabel': 'MG FxFx', 'lineColor': 'ROOT.kRed'}),
    # "Sherpa2210": sample.Sample('Sherpa2210', None, **{'add': [], 'subtract': [], 'legendLabel': 'Sherpa 2.2.10', 'lineColor': 'ROOT.kBlue'}),
}

# label dict
label_dict = {
    "meson_eta": "#eta(D)",
    "meson_pt": "p_{T}(D)",
    "lep_eta": "#eta(lep)",
}


def create_histogram(variation, x_low, x_high, y, y_up, y_dn, name):
    logging.info(f"creating histogram for {variation}")
    x_bins = x_low + [x_high[-1]]
    x_bins_array = array.array('d', x_bins)
    h = ROOT.TH1D(f"{name}_{variation}", f"{name}_{variation}", len(x_bins) - 1, x_bins_array)
    logging.info(f"created histogram {h}, filling entries now")
    integral = 0
    one_up = 0
    one_dn = 0
    for i in range(h.GetNbinsX()):
        bin_width = h.GetBinWidth(i + 1)
        integral += bin_width * y[i]
        one_up += bin_width * y_up[i]
        one_dn += bin_width * y_dn[i]
        err = (y_up[i] + y_dn[i]) / 2.
        h.SetBinContent(i + 1, y[i])
        h.SetBinError(i + 1, err)
    utils.eprint(f"variation {variation}: {integral:.3f} + {one_up:.3f} - {one_dn:.3f}")
    return h


def process_file(f, name):
    histograms = {}
    at_variation = False
    save_histo = False
    for line in f:
        if "BEGIN HISTO1D" in line and ("_VarBand" in line or "LOMG" in line or ("MGFxFx" in line and "[" not in line) or ("Sherpa2210[MUR1_MUF1_PDF30320]" in line)):  # noqa: E501
            logging.info(f"{line}")
            at_variation = True
            if "_VarBand" in line:
                variation = re.findall("BEGIN HISTO1D .([\w]+)\[_VarBand\]", line)[0]  # noqa: W605
            elif "LOMG" in line:
                variation = "LOMG"
            elif "MGFxFx" in line:
                variation = "MGFxFx"
            elif "Sherpa2210" in line:
                variation = "Sherpa2210"
            x_low = []
            x_high = []
            y = []
            y_up = []
            y_dn = []
            print(f"found variation {variation}")
        elif at_variation and "xlow" in line and "xhigh" in line and "val" in line:
            save_histo = True
        elif "END HISTO1D" in line:
            if at_variation:
                h = create_histogram(variation, x_low, x_high, y, y_up, y_dn, name)
                histograms[variation] = h
            at_variation = False
            save_histo = False
        elif save_histo:
            vals = [float(x) for x in line.split()]
            x_low += [vals[0]]
            x_high += [vals[1]]
            y += [vals[2]]
            y_up += [vals[3]]
            y_dn += [vals[4]]
    return histograms


def main(options):

    # save to file
    f = ROOT.TFile(os.path.join(options.output, "output.root"), "RECREATE")
    f.Close()

    for c in options.channels.split(","):

        # make output folder if not exist
        if not os.path.isdir(os.path.join(options.output, c)):
            os.makedirs(os.path.join(options.output, c))

        # bookkeeping for later
        first_plot = True

        for v in options.vars.split(","):
            file_name = os.path.join(options.input, f"{c}_{v}.dat")
            f = open(file_name, "r")
            logging.info(f"Succesfully opened file {file_name} as {f}")

            # check if last plot
            last_plot = v == options.vars.split(",")[-1]

            # process the dat file..
            histograms = process_file(f, f"{c}_{v}")
            print(histograms)

            # normalize to unity if requested
            if options.normalize:
                for h in histograms.values():
                    h.Scale(1. / h.GetSum())

            # var
            var = variable.Variable(v, **{"label": label_dict[v], "unit": ("GeV" if v in ["meson_pt", "lep_pt"] else "")})

            # channel
            chan = channel.Channel(c, ["W+D", c], "", [], [])

            # samples
            samples = [samples_dict[x] for x in histograms.keys() if x in samples_dict]

            # mc map
            mc_map = {samples_dict[x]: histograms[x] for x in histograms if x in samples_dict}

            # canvas
            yaxis_label = f"d#sigma / d{label_dict[v]} [pb]"
            canv = utils.make_canvas_mc_ratio(mc_map[samples[0]], var, chan, "Ratio", x=800, y=800, ratio_range=[0.31, 1.69], events=yaxis_label)

            # configure histograms
            canv.configure_histograms(mc_map, True)

            # top pad
            canv.pad1.cd()
            errors = []
            canv.pad1.cd()
            for s in samples:
                print(s.name)
                fcolor = mc_map[s].GetLineColor()
                gr_mc_stat_err, _ = utils.make_stat_err(mc_map[s])
                gr_mc_stat_err.SetLineColor(fcolor)
                gr_mc_stat_err.SetFillColorAlpha(fcolor, 0.75)
                gr_mc_stat_err.SetFillStyle(1001)
                errors += [gr_mc_stat_err]
                # if "PDF260000" in s.name or "QCD" in s.name:
                gr_mc_stat_err.Draw("e2")
                mc_map[s].Draw("hist same")

            # make legend
            canv.make_legend(mc_map, samples, print_yields=False, leg_offset=-0.15)

            # set maximum after creating legend
            canv.set_maximum([mc_map[s] for s in samples], var, mc_map[samples[0]])

            # bottom pad
            canv.pad2.cd()

            # ratio histograms
            ratios = []
            denominator = mc_map[samples[0]].Clone(f"{mc_map[samples[0]].GetName()}_denominator")
            for i in range(0, denominator.GetNbinsX() + 2):
                denominator.SetBinError(i, 0)
            for i in range(0, len(samples)):
                h = mc_map[samples[i]].Clone(f"{mc_map[samples[i]].GetName()}_ratio")
                h.Divide(denominator)
                ratios += [h]
                fcolor = mc_map[samples[i]].GetLineColor()
                gr_mc_stat_err, _ = utils.make_stat_err(h)
                gr_mc_stat_err.SetLineColor(fcolor)
                gr_mc_stat_err.SetFillColorAlpha(fcolor, 0.25)
                gr_mc_stat_err.SetFillStyle(1001)
                errors += [gr_mc_stat_err]
                # if "PDF260000" in h.GetName() or "QCD" in h.GetName():
                gr_mc_stat_err.Draw("e2")
                h.Draw("hist same")

            # save output
            h_lomg = mc_map[samples_dict["LOMG"]]
            h_nominal = mc_map[samples_dict["PDF260000"]].Clone(mc_map[samples_dict["PDF260000"]].GetName() + "_NOMINAL")
            h_pdf = mc_map[samples_dict["PDF260000"]].Clone(mc_map[samples_dict["PDF260000"]].GetName() + "_PDF")
            h_qcd = mc_map[samples_dict["QCD"]].Clone(mc_map[samples_dict["QCD"]].GetName() + "_RATIO")
            h_PDF261000 = mc_map[samples_dict["PDF261000"]].Clone(mc_map[samples_dict["PDF261000"]].GetName() + "_RATIO")
            h_PDF303600 = mc_map[samples_dict["PDF303600"]].Clone(mc_map[samples_dict["PDF303600"]].GetName() + "_RATIO")
            h_PDF25200 = mc_map[samples_dict["PDF25200"]].Clone(mc_map[samples_dict["PDF25200"]].GetName() + "_RATIO")
            h_PDF65060 = mc_map[samples_dict["PDF65060"]].Clone(mc_map[samples_dict["PDF65060"]].GetName() + "_RATIO")
            h_PDF14400 = mc_map[samples_dict["PDF14400"]].Clone(mc_map[samples_dict["PDF14400"]].GetName() + "_RATIO")
            h_PDF14600 = mc_map[samples_dict["PDF14600"]].Clone(mc_map[samples_dict["PDF14600"]].GetName() + "_RATIO")
            for i in range(0, h_nominal.GetNbinsX() + 2):
                h_pdf.SetBinContent(i, h_pdf.GetBinContent(i) + h_pdf.GetBinError(i))
                h_qcd.SetBinContent(i, h_qcd.GetBinContent(i) + h_qcd.GetBinError(i))
                h_nominal.SetBinError(i, 0)
                h_pdf.SetBinError(i, 0)
                h_qcd.SetBinError(i, 0)
                h_PDF261000.SetBinError(i, 0)
                h_PDF303600.SetBinError(i, 0)
                h_PDF25200.SetBinError(i, 0)
                h_PDF65060.SetBinError(i, 0)
                h_PDF14400.SetBinError(i, 0)
                h_PDF14600.SetBinError(i, 0)
            h_nominal.Divide(h_lomg)
            h_pdf.Divide(h_lomg)
            h_qcd.Divide(h_lomg)
            h_PDF261000.Divide(h_lomg)
            h_PDF303600.Divide(h_lomg)
            h_PDF25200.Divide(h_lomg)
            h_PDF65060.Divide(h_lomg)
            h_PDF14400.Divide(h_lomg)
            h_PDF14600.Divide(h_lomg)

            # save to file
            f = ROOT.TFile(os.path.join(options.output, "output.root"), "UPDATE")
            h_nominal.Write()
            h_pdf.Write()
            h_qcd.Write()
            h_PDF261000.Write()
            h_PDF303600.Write()
            h_PDF25200.Write()
            h_PDF65060.Write()
            h_PDF14400.Write()
            h_PDF14600.Write()
            f.Close()

            # Print out
            canv.print_all(options.output, c, v, multipage_pdf=True, first_plot=first_plot, last_plot=last_plot, as_png=False)
            first_plot = False


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-i', '--input',
                      action="store", dest="input",
                      help="analysis config file")
    parser.add_option('-o', '--output',
                      action="store", dest="output",
                      help="output folder name",
                      default="rivet_bands")
    parser.add_option('-c', '--channels',
                      action="store", dest="channels",
                      help="run over a subset of channels (comma separated)")
    parser.add_option('-v', '--vars',
                      action="store", dest="vars",
                      help="run over a subset of variables (comma separated)")
    parser.add_option('-n', '--normalize',
                      action="store_true", dest="normalize",
                      help="normalize to unity")

    # parse input arguments
    options, args = parser.parse_args()

    # output name
    out_name = options.output

    # make output folder if not exist
    if not os.path.isdir(out_name):
        os.makedirs(out_name)

    # do the plotting
    main(options)
