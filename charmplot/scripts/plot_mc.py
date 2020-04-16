#!/usr/bin/env python
from charmplot.common import utils
from charmplot.common import www
from charmplot.control import inputDataReader
import logging
import os
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


def main(options, conf, reader):

    # loop through all channels and variables
    for c in conf.channels:

        # filter channels
        if options.channels:
            if c.name not in options.channels.split(","):
                logging.debug(f"skipping channel {c.name}")
                continue

        # keep track of first/last plot of each channel
        first_plot = True

        # samples
        if not c.samples:
            logging.critical(f"no samples given for channel {c}")
            sys.exit(1)
        samples = [conf.get_sample(s) for s in c.samples]
        logging.info(f"making plots for channel {c} with samples {samples}")

        # list of variables
        variables = utils.get_variables(options, conf, reader, c, samples[0])

        # make channel folder if not exist
        if not os.path.isdir(os.path.join(options.output, c.name)):
            os.makedirs(os.path.join(options.output, c.name))
        for v in variables:

            # check if last plot
            last_plot = v == variables[-1]

            # mc map
            mc_map = {s: reader.get_histogram(s, c, conf.get_var(v)) for s in samples}

            # canvas
            canv = utils.make_canvas(mc_map[samples[0]], conf.get_var(v), c, x=800, y=800, y_split = 0)

            # configure histograms
            canv.configure_histograms(mc_map)

            # top pad
            canv.pad1.cd()
            for s in samples:
                mc_map[s].Draw("hist same")

            # make legend
            canv.make_legend(data = None, mc_map = mc_map, samples = samples)

            # set maximum after creating legend
            canv.set_maximum([mc_map[s] for s in samples], conf.get_var(v), mc_map[samples[0]])

            # Print out
            canv.print_all(options.output, c.name, v, multipage_pdf=True, first_plot=first_plot, last_plot=last_plot, as_png=options.stage_out)
            first_plot = False


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-a', '--analysis-config',
                      action="store", dest="analysis_config",
                      help="analysis config file")
    parser.add_option('-c', '--channels',
                      action="store", dest="channels",
                      help="run over a subset of channels (comma separated)")
    parser.add_option('-v', '--vars',
                      action="store", dest="vars",
                      help="run over a subset of variables (comma separated)")
    parser.add_option('-n', '--normalize',
                      action="store_true", dest="normalize",
                      help="normalize to luminosity")
    parser.add_option('-o', '--output-file',
                      action="store", dest="output",
                      help="save histograms to an output file")
    parser.add_option('--stage-out',
                      action="store_true", dest="stage_out",
                      help="copy plots to the www folder")

    # parse input arguments
    options, args = parser.parse_args()

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

    # stage-out to the www folder
    if options.stage_out:
        www.stage_out_plots(options.output, conf.get_variables(), x=300, y=300)
