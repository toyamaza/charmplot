# flake8: noqa: C901
#!/usr/bin/env python
from charmplot.common import utils
from charmplot.common import www
from charmplot.control import globalConfig
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


def read_trex_input(trex_folder):
    trex_histograms = {}
    for subdir, dirs, files in os.walk(os.path.join(trex_folder, "Histograms")):
        for file in files:
            if "postFit" in file:
                channel = file.replace("_postFit.root", "")
                trex_histograms[channel] = os.path.join(subdir, file)
    logging.info(f"Found TREx histograms {trex_histograms}")
    return trex_histograms


def make_trex_folder(trex_folder):
    if not os.path.isdir(trex_folder):
        os.makedirs(trex_folder)


def mass_fit(conf, reader, c, samples):
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


def main(options, conf, reader):

    # make the output file only once
    histogram_file_made = False

    # trex histogram files (one per sample)
    trex_histograms = {}

    # read trex input
    trex_post_fit_histograms = {}
    if options.trex_input:
        trex_post_fit_histograms = read_trex_input(options.trex_input)

    # loop through all channels and variables
    for c in conf.channels:

        # output root file
        if c.save_to_file and not histogram_file_made:
            histogram_file_made = True
            out_file_name = os.path.join(conf.out_name, "histograms.root")
            out_file = ROOT.TFile(out_file_name, "RECREATE")
            out_file.Close()

        # filter channels
        if options.channels:
            if c.name not in options.channels.split(","):
                logging.debug(f"skipping channel {c.name}")
                continue

        # skip channels
        if not c.make_plots and not c.save_to_file:
            continue

        # make channel folder if not exist
        if not os.path.isdir(os.path.join(conf.out_name, c.name)):
            os.makedirs(os.path.join(conf.out_name, c.name))

        # used MC samples in channel (default or channel specific)
        samples = utils.get_samples(conf, c)

        # trex output
        if options.trex:
            trex_folder = os.path.join(conf.out_name, options.trex)
            for s in samples + [conf.get_data()]:
                make_trex_folder(trex_folder)
                sample_file_name = os.path.join(trex_folder, f"{s.shortName}.root")
                if s.shortName not in trex_histograms:
                    sample_out_file = ROOT.TFile(sample_file_name, "RECREATE")
                    trex_histograms[s.shortName] = sample_file_name
                    sample_out_file.Close()

        # perform likelihood fit
        fit = utils.likelihood_fit(conf, reader, c, samples)

        # keep track of first/last plot of each channel
        first_plot = True

        # list of variables
        variables = utils.get_variables(options, conf, reader, c)

        # mass fit
        if c.mass_fit:
            mass_fit(conf, reader, c, samples)

        for v in variables:

            # variable object
            var = conf.get_var(v)

            # systematics
            systematics = [None]
            for sys in systematics:

                # check if last plot
                last_plot = v == variables[-1] and sys == systematics[-1]

                # data histogram
                h_data = reader.get_histogram(conf.get_data(), c, var)
                if not h_data:
                    continue

                # read input MC histograms (and scale them)
                mc_map = utils.read_samples(conf, reader, c, var, samples, fit, force_positive=c.force_positive, sys=sys)
                if not mc_map:
                    continue

                # trex post-fit
                trex_mc_tot = None
                trex_mc_stat_err = None
                trex_mc_stat_err_only = None
                if trex_post_fit_histograms and c.trex_subtraction:
                    channelOS = conf.get_channel(c.trex_subtraction["OS"])
                    channelSS = conf.get_channel(c.trex_subtraction["SS"])
                    h_data, trex_mc_tot, trex_mc_stat_err, trex_mc_stat_err_only = utils.trex_subtraction(channelOS, channelSS, var, mc_map, trex_post_fit_histograms)
                    c.label += ["post-fit"]
                elif c.name in trex_post_fit_histograms:
                    h_data, trex_mc_tot, trex_mc_stat_err, trex_mc_stat_err_only = utils.read_trex_input(c, var, mc_map, trex_post_fit_histograms)
                    c.label += ["post-fit"]

                # continue if not make plots
                if not c.make_plots or not var.make_plots:
                    continue

                # scale factors for this channel (only for dispaly)
                scale_factors = utils.read_scale_factors(c.scale_factors)

                # canvas
                canv = utils.make_canvas(h_data, var, c, x=800, y=800, fit=fit, scale_factors=scale_factors, sys=sys)

                # configure histograms
                canv.configure_histograms(mc_map, h_data)

                # save histograms to root file
                if c.save_to_file:
                    utils.save_to_file(out_file_name, c, var, h_data, mc_map)
                    if options.trex:
                        utils.save_to_trex_file(trex_folder, c, var, h_data, mc_map, trex_histograms, sys=sys)

                # stack and total mc
                hs = utils.make_stack(samples, mc_map)
                if trex_mc_tot:
                    h_mc_tot = trex_mc_tot
                else:
                    h_mc_tot = utils.make_mc_tot(hs, f"{c}_{v}_mc_tot")

                # ratio
                h_ratio = utils.make_ratio(h_data, h_mc_tot)

                # mc error
                if trex_mc_stat_err and trex_mc_stat_err_only:
                    gr_mc_stat_err, gr_mc_stat_err_only = trex_mc_stat_err, trex_mc_stat_err_only
                    if c.name not in ["OS-SS_2018_el_SR_2tag_Dplus", "OS-SS_2018_mu_SR_2tag_Dplus"]:
                        canv.set_ratio_range(0.89, 1.11, override=True)
                else:
                    gr_mc_stat_err, gr_mc_stat_err_only = utils.make_stat_err(h_mc_tot)

                # top pad
                canv.pad1.cd()
                hs.Draw("same hist")
                h_mc_tot.Draw("same hist")
                gr_mc_stat_err.Draw("e2")
                h_data.Draw("same pe")

                # make legend
                canv.make_legend(h_data, h_mc_tot, mc_map, samples, print_yields=True, show_error=(trex_mc_tot is None))

                # set maximum after creating legend
                canv.set_maximum((h_data, h_mc_tot), var, mc_min=utils.get_mc_min(mc_map, samples))

                # bottom pad
                canv.pad2.cd()
                if not c.qcd_template:
                    gr_mc_stat_err_only.Draw("le2")
                    h_ratio.Draw("same pe")
                else:
                    h_qcd_frac, h_qcd_frac_err = utils.get_fraction_histogram(mc_map[conf.get_sample(c.qcd_template)], h_data)
                    canv.draw_qcd_frac(h_qcd_frac, h_qcd_frac_err)

                # Print out
                canv.print_all(conf.out_name, c.name, v, multipage_pdf=True, first_plot=first_plot, last_plot=last_plot, as_png=options.stage_out)
                first_plot = False

                # close output file
                if c.save_to_file:
                    out_file.Close()


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
    parser.add_option('--suffix',
                      action="store", dest="suffix",
                      help="suffix for the output name")
    parser.add_option('--stage-out',
                      action="store_true", dest="stage_out",
                      help="copy plots to the www folder")
    parser.add_option('--trex',
                      action="store", dest="trex",
                      help="make output histograms for TRExFitter")
    parser.add_option('--trex-input',
                      action="store", dest="trex_input",
                      help="import post-fit trex plots")

    # parse input arguments
    options, args = parser.parse_args()

    # analysis configs
    config = options.analysis_config

    # output name
    out_name = config
    if options.suffix:
        out_name = out_name.split("/")
        out_name[0] += "_" + options.suffix
        out_name = "/".join(out_name)

    # make output folder if not exist
    if not os.path.isdir(out_name):
        os.makedirs(out_name)

    # read inputs
    conf = globalConfig.GlobalConfig(config, out_name)
    reader = inputDataReader.InputDataReader(conf)

    # do the plotting
    main(options, conf, reader)

    # stage-out to the www folder
    if options.stage_out:
        www.stage_out_plots(conf.out_name, conf.get_variables(), x=300, y=300)
