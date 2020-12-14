# flake8: noqa: C901
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

    # make the output file only once
    histogram_file_made = False

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

        # keep track of first/last plot of each channel
        first_plot = True

        # skip channels
        if not (c.make_plots or c.save_to_file):
            continue

        # samples
        if not c.samples:
            logging.critical(f"no samples given for channel {c}")
            sys.exit(1)
        samples = [conf.get_sample(s) for s in c.samples]
        logging.info(f"making plots for channel {c} with samples {samples}")

        # list of variables
        variables = utils.get_variables(options, conf, reader, c, samples[0])

        # make channel folder if not exist
        if not os.path.isdir(os.path.join(options.analysis_config, c.name)):
            os.makedirs(os.path.join(options.analysis_config, c.name))
        for v in variables:

            # variable object
            var = conf.get_var(v)

            # check if last plot
            last_plot = v == variables[-1]

            # read input MC histograms (and scale them)
            mc_map = utils.read_samples(conf, reader, c, var, samples)
            if not mc_map:
                continue

            # save histograms to root file
            if c.save_to_file:
                utils.save_to_file(out_file_name, c, var, None, mc_map)

            # continue if not make plots
            if not c.make_plots or not var.make_plots:
                continue

            # canvas
            canv = utils.make_canvas(mc_map[samples[0]], var, c, x=800, y=800, y_split=0)

            # configure histograms
            canv.configure_histograms(mc_map)

            # top pad
            hists = [mc_map[s] for s in mc_map]
            errors = []
            canv.pad1.cd()
            if not options.do_stack:
                for s in samples:
                    if s not in mc_map.keys():
                        continue
                    fcolor = mc_map[s].GetLineColor()
                    gr_mc_stat_err, gr_mc_stat_err_only = utils.make_stat_err(mc_map[s])
                    gr_mc_stat_err.SetLineColor(fcolor)
                    gr_mc_stat_err.SetFillColorAlpha(fcolor, 0.25)
                    gr_mc_stat_err.SetFillStyle(1001)
                    errors += [gr_mc_stat_err]
                    gr_mc_stat_err.Draw("e2")
                    mc_map[s].Draw("hist same")
            else:
                hs = utils.make_stack(samples, mc_map)
                hs.Draw("samehist")

            # make legend
            canv.make_legend(data=None, mc_map=mc_map, samples=samples, print_yields=True)

            # set maximum after creating legend
            if not options.do_stack:
                canv.set_maximum(hists, var, utils.get_mc_min(mc_map, samples))
            else:
                h_mc_tot = utils.make_mc_tot(hs, "mc_tot")
                canv.set_maximum([h_mc_tot], var, mc_min=utils.get_mc_min(mc_map, samples))

            # Print out
            canv.print_all(options.analysis_config, c.name, v, multipage_pdf=True, first_plot=first_plot, last_plot=last_plot, as_png=options.stage_out)
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
    parser.add_option('--suffix',
                      action="store", dest="suffix",
                      help="suffix for the output name")
    parser.add_option('--stage-out',
                      action="store_true", dest="stage_out",
                      help="copy plots to the www folder")
    parser.add_option('--stack',
                      action="store_true", dest="do_stack",
                      help="make stacked plots of categories")

    # parse input arguments
    options, args = parser.parse_args()

    # output name
    out_name = options.analysis_config
    if options.suffix:
        out_name = out_name.split("/")
        out_name[0] += "_" + options.suffix
        out_name = "/".join(out_name)

    # make output folder if not exist
    if not os.path.isdir(out_name):
        os.makedirs(out_name)

    # read inputs
    from charmplot.control import globalConfig
    conf = globalConfig.GlobalConfig(options.analysis_config, out_name)
    reader = inputDataReader.InputDataReader(conf)

    # do the plotting
    main(options, conf, reader)

    # stage-out to the www folder
    if options.stage_out:
        www.stage_out_plots(options.analysis_config, conf.get_variables(), x=300, y=300)
