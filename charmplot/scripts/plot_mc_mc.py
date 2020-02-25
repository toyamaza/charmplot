#!/usr/bin/env python
from charmplot.common import utils
from charmplot.control import inputDataReader
import logging
import os
import ROOT
import sys

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasLabels.C"))
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasUtils.C"))
ROOT.SetAtlasStyle()

# logging
root = logging.getLogger()
root.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
root.addHandler(handler)


def main(options, conf, reader):

    # loop through all channels and variables
    for c in conf.channels:

        # keep track of first/last plot of each channel
        first_plot = True

        # list of variables
        variables = utils.get_variables(options, conf, reader, c)

        # make channel folder if not exist
        if not os.path.isdir(os.path.join(options.output, c.name)):
            os.makedirs(os.path.join(options.output, c.name))
        for v in variables:

            # check if last plot
            last_plot = v == variables[-1]

            # mc samples
            samples = [conf.get_sample(s) for s in options.samples.split(",")]
            if not samples or len(samples) % 2 != 0:
                logging.critical(f"incorect input for samples: {samples}")
                sys.exit(1)

            # mc map
            mc_map = {s: reader.get_histogram(s, c, v) for s in samples}

            # canvas
            canv = utils.make_canvas_mc_ratio(mc_map[samples[0]], conf.get_var(v), c)

            # configure histograms
            canv.configure_histograms(mc_map, conf.get_var(v))

            # top pad
            canv.pad1.cd()
            for s in samples:
                mc_map[s].Draw("hist same")

            # make legend
            canv.make_legend(mc_map, samples)

            # set maximum after creating legend
            canv.set_maximum([mc_map[s] for s in samples], conf.get_var(v))

            # bottom pad
            canv.pad2.cd()
            canv.set_ratio_range(0, 1.99)

            # ratio histograms
            ratios = []
            for i in range(0, len(samples), 2):
                h = mc_map[samples[i]].Clone(f"{mc_map[samples[i]].GetName()}_ratio")
                h.Divide(mc_map[samples[i + 1]])
                ratios += [h]
                h.Draw("hist same")

            # Print out
            canv.print_all(options.output, c.name, v, multipage_pdf=True, first_plot=first_plot, last_plot=last_plot)
            first_plot = False


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-c', '--analysis-config',
                      action="store", dest="analysis_config",
                      help="analysis config file")
    parser.add_option('-v', '--vars',
                      action="store", dest="vars",
                      help="run over a subset of variables (comma separated)")
    parser.add_option('-o', '--output-file',
                      action="store", dest="output",
                      help="save histograms to an output file")
    parser.add_option('-s', '--samples',
                      action="store", dest="samples",
                      help="comma separated list of samples to compare")

    # parse input arguments
    options, args = parser.parse_args()

    # mandatory
    assert options.samples

    # output file
    if not options.output:
        options.output = options.analysis_config

    # make output folder if not exist
    if not os.path.isdir(options.output):
        os.makedirs(options.output)

    # read inputs
    from charmplot.control import globalConfig
    conf = globalConfig.GlobalConfig(options.analysis_config)
    reader = inputDataReader.InputDataReader(conf)

    # do the plotting
    main(options, conf, reader)
