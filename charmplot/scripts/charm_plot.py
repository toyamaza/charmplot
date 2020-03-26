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

        # skip channels
        if not c.make_plots:
            continue

        # filter channels
        if options.channels:
            if c.name not in options.channels.split(","):
                logging.debug(f"skipping channel {c.name}")
                continue

        # used MC samples in channel (default or channel specific)
        samples = utils.get_samples(conf, c)

        # perform likelihood fit
        fit = utils.likelihood_fit(conf, reader, c, samples)

        # scale factors for this channel
        scale_factors = utils.read_scale_factors(options.output, c.scale_factors)

        # keep track of first/last plot of each channel
        first_plot = True

        # list of variables
        variables = utils.get_variables(options, conf, reader, c)

        # mass fit
        if c.mass_fit:
            fit_var = conf.get_var(c.mass_fit["var"])
            h_data = reader.get_histogram(conf.get_data(), c, fit_var)
            mc_map = {}
            for s in samples:
                mc_map[s] = reader.get_histogram(s, c, fit_var)
            hs = utils.make_stack(samples, mc_map)
            h_mc_tot = utils.make_mc_tot(hs, f"{c}_{fit_var.name}_mc_tot")
            # fit data
            utils.mass_fit(conf, h_data, c, "Data")
            # fit MC
            utils.mass_fit(conf, h_mc_tot, c, "MC")

        # make channel folder if not exist
        if not os.path.isdir(os.path.join(options.output, c.name)):
            os.makedirs(os.path.join(options.output, c.name))
        for v in variables:

            # check if last plot
            last_plot = v == variables[-1]

            # data histogram
            h_data = reader.get_histogram(conf.get_data(), c, conf.get_var(v))

            # read input MC histograms (and scale them)
            mc_map = {}
            for s in samples:
                # read MC histogram
                h = reader.get_histogram(s, c, conf.get_var(v))
                mc_map[s] = h

                # scale histogram if performed likelihood fit
                if fit:
                    h.Scale(fit.result[s.name][0])

                # scale histogram if given input scale factors
                # TODO: improve the sample name propagation
                if scale_factors:
                    if s.name in scale_factors:
                        sf = scale_factors[s.name]
                    else:
                        sf = scale_factors[c.scale_factors['scale_factors'][s.name]]
                    h.Scale(sf[0])

            # canvas
            canv = utils.make_canvas(h_data, conf.get_var(v), c, x=800, y=800, fit=fit)

            # configure histograms
            canv.configure_histograms(mc_map, h_data, conf.get_var(v))

            # stack and total mc
            hs = utils.make_stack(samples, mc_map)
            h_mc_tot = utils.make_mc_tot(hs, f"{c}_{v}_mc_tot")

            # ratio
            h_ratio = utils.make_ratio(h_data, h_mc_tot)

            # mc error
            gr_mc_stat_err, gr_mc_stat_err_only = utils.make_stat_err(h_mc_tot)

            # top pad
            canv.pad1.cd()
            hs.Draw("same hist")
            h_mc_tot.Draw("same hist")
            gr_mc_stat_err.Draw("e2")
            h_data.Draw("same pe")

            # make legend
            canv.make_legend(h_data, h_mc_tot, mc_map, samples, print_yields=True)

            # set maximum after creating legend
            canv.set_maximum((h_data, h_mc_tot), conf.get_var(v), mc_min=mc_map[samples[-1]])

            # bottom pad
            canv.pad2.cd()
            if not c.qcd_template:
                gr_mc_stat_err_only.Draw("le2")
                h_ratio.Draw("same pe")
            else:
                h_qcd_frac, h_qcd_frac_err = utils.get_fraction_histogram(mc_map[conf.get_sample(c.qcd_template)], h_data)
                canv.draw_qcd_frac(h_qcd_frac, h_qcd_frac_err)

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
