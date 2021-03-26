#!/usr/bin/env python
from charmplot.common import utils
from charmplot.common import www
from charmplot.control import globalConfig
from charmplot.control import inputDataReader
from multiprocessing import Pool
import logging
import os
import ROOT
import sys
import time

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


def wrapper(args,):
    options, conf, c = args
    process_channel(options, conf, c)


def process_channel(options, conf, c):

    # continue if proxy channel
    if not (c.save_to_file or c.make_plots):
        return

    logging.info(f"start processing channel {c.name}")

    # read inputs
    reader = inputDataReader.InputDataReader(conf)

    # read trex input
    trex_post_fit_histograms = {}
    if options.trex_input:
        trex_post_fit_histograms = read_trex_input(options.trex_input)

    # trex histogram files (one per sample)
    trex_histograms = {}
    if (options.trex):
        trex_folder = os.path.join(conf.out_name, options.trex)

    # output root file
    if c.save_to_file:
        out_file_name = os.path.join(conf.out_name, f"histograms_tmp_{c.name}.root")
        out_file = ROOT.TFile(out_file_name, "RECREATE")
        out_file.Close()

    # make channel folder if not exist
    if not os.path.isdir(os.path.join(conf.out_name, c.name)):
        os.makedirs(os.path.join(conf.out_name, c.name))

    # used MC samples in channel (default or channel specific)
    samples = utils.get_samples(conf, c)

    # trex output
    if options.trex and c.make_plots:
        for s in samples + [conf.get_data()]:
            sample_file_name = os.path.join(trex_folder, f"{s.shortName}_tmp_{c.name}.root")
            logging.info(f"created file for trex input {sample_file_name}")
            if s.shortName not in trex_histograms:
                sample_out_file = ROOT.TFile(sample_file_name, "RECREATE")
                trex_histograms[s.shortName] = sample_file_name
                sample_out_file.Close()
                if s.makeGhostSample:
                    sample_file_name_ghost = os.path.join(trex_folder, f"{s.shortName}_GHOST_tmp_{c.name}.root")
                    sample_out_file = ROOT.TFile(sample_file_name_ghost, "RECREATE")
                    trex_histograms[s.shortName + "_GHOST"] = sample_file_name_ghost
                    sample_out_file.Close()

    # perform likelihood fit
    fit = utils.likelihood_fit(conf, reader, c, samples)

    # keep track of first/last plot of each channel
    first_plot = True

    # list of variables
    variables = utils.get_variables(options, conf, reader, c)
    assert len(variables) > 0

    # systematics
    systematics = conf.get_systematics()

    # no sys if TREx import
    if options.trex_input:
        systematics = None

    for v in variables:

        # variable object
        var = conf.get_var(v)

        # check if last plot
        last_plot = v == variables[-1]

        # data histogram
        h_data = reader.get_histogram(conf.get_data(), c, var)
        if not h_data:
            continue

        # read input MC histograms (and scale them)
        mc_map = utils.read_samples(conf, reader, c, var, samples, fit, force_positive=c.force_positive)
        if not mc_map:
            continue

        # systematics histograms
        if systematics:
            mc_map_sys = utils.read_sys_histograms(conf, reader, c, var, samples, fit, systematics, mc_map)

        # trex post-fit
        trex_mc_tot = None
        trex_mc_stat_err = None
        trex_mc_stat_err_only = None
        if trex_post_fit_histograms and c.trex_subtraction:
            channelOS = conf.get_channel(c.trex_subtraction["OS"])
            channelSS = conf.get_channel(c.trex_subtraction["SS"])
            h_data, trex_mc_tot, trex_mc_stat_err, trex_mc_stat_err_only = utils.trex_subtraction(
                channelOS, channelSS, var, mc_map, trex_post_fit_histograms)
            c.label += ["post-fit"]
        elif c.name in trex_post_fit_histograms:
            h_data, trex_mc_tot, trex_mc_stat_err, trex_mc_stat_err_only = utils.read_trex_input(c, var, mc_map, trex_post_fit_histograms)
            c.label += ["post-fit"]

        # replace samples
        for replace_sample, replace_channel in c.replacement_samples.items():
            logging.info(f"replacing sample {replace_sample} with channel {replace_channel}")
            utils.replace_sample(conf, mc_map, reader, c, var, replace_sample, replace_channel, mc_map_sys if systematics else None)

        # save histograms to root file
        if c.save_to_file:
            utils.save_to_file(out_file_name, c, var, h_data, mc_map)

        # continue if not make plots
        if not c.make_plots or not var.make_plots:
            continue

        # canvas
        canv = utils.make_canvas(h_data, var, c, x=800, y=800, fit=fit)

        # configure histograms
        canv.configure_histograms(mc_map, h_data, style=conf.style)

        # save histograms to root file
        if options.trex:
            utils.save_to_trex_file(trex_folder, c, var, h_data, mc_map, trex_histograms)
            if systematics:
                for group in systematics:
                    variations = systematics[group].get('variations')
                    affecting = systematics[group].get('affecting')
                    sys_type = systematics[group].get('type')
                    if sys_type in ['updown', 'alt_sample', 'overall']:
                        for syst in variations:
                            utils.save_to_trex_file(trex_folder, c, var, None, mc_map_sys[group][syst], trex_histograms, syst, affecting)
                    elif sys_type == 'minmax':
                        for s in mc_map:
                            if affecting and s.shortName not in affecting:
                                continue
                            sys_hists = [mc_map_sys[group][syst][s] for syst in variations]
                            nominal_hist = mc_map[s]
                            gr_mc_sys_err, _ = utils.make_minmax_err(nominal_hist, sys_hists)
                            h_up, _ = utils.sys_graph_to_hists(gr_mc_sys_err, nominal_hist, group)
                            utils.save_histograms_to_trex_file(trex_folder, c, var, h_up, s, trex_histograms, f"{group}_up")
                    else:
                        for s in mc_map:
                            if affecting and s.shortName not in affecting:
                                continue
                            sys_hists = [mc_map_sys[group][syst][s] for syst in variations]
                            nominal_hist = mc_map[s]
                            gr_mc_sys_err, _ = utils.make_pdf_err(nominal_hist, sys_hists, sys_type)
                            h_up, h_dn = utils.sys_graph_to_hists(gr_mc_sys_err, nominal_hist, group)
                            utils.save_histograms_to_trex_file(trex_folder, c, var, h_up, s, trex_histograms, f"{group}_up")
                            utils.save_histograms_to_trex_file(trex_folder, c, var, h_dn, s, trex_histograms, f"{group}_dn")

        # stack and total mc
        hs = utils.make_stack(samples, mc_map)
        if trex_mc_tot:
            h_mc_tot = trex_mc_tot
        else:
            h_mc_tot = utils.make_mc_tot(hs, f"{c}_{v}_mc_tot")

        # MC tot for systematics
        if systematics:
            h_mc_tot_sys = {}
            for group in systematics:
                variations = systematics[group]['variations']
                h_mc_tot_sys[group] = [utils.make_mc_tot(utils.make_stack(samples, mc_map_sys[group][syst]),
                                                         f"{c}_{v}_{group}_{syst}_mc_tot") for syst in variations]

        # ratio
        h_ratio = utils.make_ratio(h_data, h_mc_tot)

        # mc stat error
        if trex_mc_stat_err and trex_mc_stat_err_only:
            gr_mc_stat_err, gr_mc_stat_err_only = trex_mc_stat_err, trex_mc_stat_err_only
        else:
            gr_mc_stat_err, gr_mc_stat_err_only = utils.make_stat_err(h_mc_tot)

        # systematics error bands
        gr_mc_sys_err_map = []
        gr_mc_sys_err_only_map = []
        if systematics:
            for group in systematics:
                sys_type = systematics[group]['type']
                if sys_type in ['updown', 'alt_sample', 'overall']:
                    gr_mc_sys_err, gr_mc_sys_err_only = utils.make_sys_err(h_mc_tot, h_mc_tot_sys[group])
                elif sys_type == 'minmax':
                    gr_mc_sys_err, gr_mc_sys_err_only = utils.make_minmax_err(h_mc_tot, h_mc_tot_sys[group])
                else:
                    print(sys_type, " ", len(h_mc_tot_sys[group]))
                    gr_mc_sys_err, gr_mc_sys_err_only = utils.make_pdf_err(h_mc_tot, h_mc_tot_sys[group], sys_type)
                gr_mc_sys_err_map += [gr_mc_sys_err]
                gr_mc_sys_err_only_map += [gr_mc_sys_err_only]

        # total error
        gr_mc_tot_err = utils.combine_error_multiple([gr_mc_stat_err] + gr_mc_sys_err_map)
        gr_mc_tot_err_only = utils.combine_error_multiple([gr_mc_stat_err_only] + gr_mc_sys_err_only_map)

        # top pad
        canv.pad1.cd()
        hs.Draw("same hist")
        h_mc_tot.Draw("same hist")
        gr_mc_tot_err.Draw("e2")
        # gr_mc_stat_err.Draw("e2")
        h_data.Draw("same pe")

        # make legend
        canv.make_legend(h_data, h_mc_tot, mc_map, samples, print_yields=True, show_error=False)

        # set maximum after creating legend
        canv.set_maximum((h_data, h_mc_tot), var, mc_min=utils.get_mc_min(mc_map, samples))

        # find minimum
        if options.nology:
            min_negative = {}
            for s in samples:
                if s not in mc_map:
                    continue
                h = mc_map[s]
                for i in range(1, h.GetNbinsX() + 1):
                    if h.GetBinContent(i) < 0:
                        if i not in min_negative:
                            min_negative[i] = 0
                        min_negative[i] += h.GetBinContent(i)
            if min_negative.values():
                min_negative = min(min_negative.values())
                canv.proxy_up.SetMinimum(min_negative)

        # bottom pad
        canv.pad2.cd()
        if not c.qcd_template:
            gr_mc_tot_err_only.Draw("le2")
            # gr_mc_stat_err_only.Draw("le2")
            h_ratio.Draw("same pe")
        else:
            h_qcd_frac, h_qcd_frac_err = utils.get_fraction_histogram(mc_map[conf.get_sample(c.qcd_template)], h_data)
            canv.draw_qcd_frac(h_qcd_frac, h_qcd_frac_err)

        # Print out
        canv.print_all(conf.out_name, c.name, v, multipage_pdf=True, first_plot=first_plot,
                       last_plot=last_plot, as_png=options.stage_out, logy=(not options.nology))
        first_plot = False

        # close output file
        if c.save_to_file:
            out_file.Close()

    logging.info(f"finished processing channel {c.name}")


def main(options, conf):

    # trex output
    if options.trex:
        trex_folder = os.path.join(conf.out_name, options.trex)
        make_trex_folder(trex_folder)

    # make output folder
    if not os.path.isdir(conf.out_name):
        os.makedirs(conf.out_name)

    # loop through all channels and variables
    concurrnet_jobs = []
    for c in conf.channels:

        # filter channels
        if options.channels:
            if c.name not in options.channels.split(","):
                logging.debug(f"skipping channel {c.name}")
                continue

        # skip channels
        if not c.make_plots and not c.save_to_file:
            logging.debug(f"skipping channel {c.name}")
            continue

        concurrnet_jobs += [[options, conf, c]]

    print(len(concurrnet_jobs))

    # multiprocessing pool
    p = Pool(int(options.threads))
    for i, _ in enumerate(p.imap_unordered(wrapper, concurrnet_jobs)):
        print("done processing job %s/%s" % (i + 1, len(concurrnet_jobs)))


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
    parser.add_option('--nology',
                      action="store_true", dest="nology",
                      help="no log-y plots")
    parser.add_option('-t', '--threads',
                      action="store", dest="threads",
                      help="number of threads",
                      default=8)

    # parse input arguments
    options, args = parser.parse_args()

    # analysis configs
    config = options.analysis_config

    # output name
    out_name = config.replace(".yaml", "").replace(".yml", "")
    if options.suffix:
        out_name = out_name.split("/")
        if out_name[0] != ".":
            out_name[0] += "_" + options.suffix
            out_name = "/".join(out_name)
        else:
            out_name[1] += "_" + options.suffix
            out_name = "/".join(out_name)

    # config object
    conf = globalConfig.GlobalConfig(config, out_name)

    # make output folder if not exist
    if not os.path.isdir(out_name):
        os.makedirs(out_name)

    # do the plotting
    main(options, conf)

    # wait
    time.sleep(1)

    # merge trex output
    if options.trex:
        samples = {}
        trex_folder = os.path.join(conf.out_name, options.trex)
        for r, d, f in os.walk(trex_folder):
            for file in f:
                if '.root' in file and '_tmp_' in file:
                    sample = file.split("_tmp_")[0]
                    if sample not in samples:
                        samples[sample] = []
                    samples[sample] += [os.path.join(trex_folder, file)]
        jobs = [[os.path.join(trex_folder, f"{sample}.root")] + samples[sample] for sample in samples]
        print(jobs)
        p = Pool(1)
        for i, _ in enumerate(p.imap_unordered(utils.hadd_wrapper, jobs)):
            print("done processing hadd job %s/%s" % (i + 1, len(jobs)))
        os.system(f"rm {trex_folder}/*_tmp_*")

    # merge regular output
    out_files = []
    out_folder = os.path.join(conf.out_name)
    files_exist = False
    for r, d, f in os.walk(out_folder):
        for file in f:
            if '.root' in file and '_tmp_' in file:
                out_files += [os.path.join(out_folder, file)]
                files_exist = True
    if files_exist:
        jobs = [[os.path.join(out_folder, f"histograms.root")] + out_files]
        p = Pool(1)
        for i, _ in enumerate(p.imap_unordered(utils.hadd_wrapper, jobs)):
            print("done processing hadd job %s/%s" % (i + 1, len(jobs)))
        os.system(f"rm {out_folder}/*_tmp_*")

    # stage-out to the www folder
    if options.stage_out:
        www.stage_out_plots(conf.out_name, conf.get_variables(), x=300, y=300)
