#!/usr/bin/env python
import gen.utils.templates as templates
import sys
import yaml


def main(options):

    for sample_config in options.samples_config.split(","):

        # TODO: make configurable
        extra_rebin = float(options.extra_rebin)

        # TODO: make configurable?
        make_os_ss = True
        make_os_minus_ss = not options.fit_only
        os_only = options.fit_only
        force_positive = False
        # force_positive = options.fit_only or (options.fit_type == "OS/SS")

        # sample type
        if options.samples.lower() == 'truth':
            if options.fit_only:
                samples = templates.WDTruthSamplesNew(os_minus_ss_fit_configuration=(options.fit_type == "OS-SS"),
                                                      OS_and_SS_fit=(options.fit_type == "OS/SS"), MockMC=True, decayMode=options.decay_mode)
            else:
                samples = templates.WDTruthSamplesNew(OS_and_SS_fit=(options.fit_type == "OS/SS"), MockMC=False, decayMode=options.decay_mode)
            if options.replacement_samples:
                samples_replacement = templates.ReplacementSamples()
        elif options.samples.lower() == 'truth_old':
            samples = templates.WDTruthSamples()
        elif options.samples.lower() == 'flavor' or options.samples.lower() == 'flavour':
            samples = templates.WDFlavourSamples()
        elif options.samples.lower() == 'fit':
            samples = templates.WDFitSamples()
        elif options.samples.lower() == 'spg_comparison':
            samples = templates.SPGComparison()
        elif options.samples.lower() == 'bkg_comparison':
            samples = templates.BKGComparison()
        elif options.samples.lower() == 'wplusd_comparison':
            samples = templates.WDComparisonSamples()
        elif options.samples.lower() == 'flavor_comparison':
            samples = templates.WDFlavourComparison()
        elif options.samples.lower() == 'background_comparison':
            samples = templates.WDBackgroundComparison()
        elif options.samples.lower() == 'multijet_comparison':
            samples = templates.MultiJetComparison()
        elif options.samples.lower() == 'multijet_composition':
            samples = templates.MultiJetComposition()
        else:
            print(f"ERROR: Unknown samples type {options.samples}")
            sys.exit(1)

        # TODO: make configurable?
        signs = ['OS', 'SS']
        lumi = ['2015', '2016A-B', '2016C-L', '2017', '2018']
        years = []  # to be used without 'split_by_period' option in charmpp
        leptons = ['el', 'mu']
        charges = ['minus', 'plus']
        btags = ['0tag', '1tag']
        ptbins = ['']
        # ptbins = ['pt_bin1', 'pt_bin2', 'pt_bin3', 'pt_bin4', 'pt_bin5', '']

        # binning for the 1-tag region
        # use '40' for 1 bin for D+
        btag_bin = 2

        # replace samples to increase mc stats
        if options.replacement_samples:
            replacement_samples = {
                'Wjets_emu_Matched': '<charge>_SPG_Matched',
                'Wjets_emu_411MisMatched': '<charge>_SPG_411MisMatched',
                'Wjets_emu_Charm': '<charge>_SPG_CharmMisMatched',
                'Wjets_emu_Rest': '<charge>_SPG_Wjets_emu_Rest',
                'Wjets_emu_MisMatched': '<charge>_SPG_Wjets_emu_MisMatched',
            }
        else:
            replacement_samples = {}

        # systematics
        systematics = []
        if options.systematics:
            systematics = [
                'experimental',
                'matrix_method',
                'proxy_norm',
                'ttbar_theory_alt_samples',
                'ttbar_theory_choice',
                'ttbar_theory_pdf',
                'ttbar_theory_qcd',
                'wjets_theory',
            ]

        # Base config
        config = templates.DataMCConfig(variables=options.variables,
                                        sample_config=sample_config,
                                        systematics=systematics).to_dict()

        # Helper object to generate channels
        channelGenerator = templates.ChannelGenerator(config=config,
                                                      samples=samples,
                                                      decay_mode=options.decay_mode,
                                                      signs=signs,
                                                      years=years,
                                                      leptons=leptons,
                                                      charges=charges,
                                                      btags=btags,
                                                      ptbins=ptbins,
                                                      process_string=options.process_string,
                                                      sample_config=sample_config,
                                                      force_positive=force_positive,
                                                      replacement_samples=replacement_samples,
                                                      os_minus_ss_fit_configuration=(options.fit_type == "OS-SS"))

        # Helper object to generate channels
        if options.replacement_samples:
            channelGeneratorReplacement = templates.ChannelGenerator(config=config,
                                                                     samples=samples_replacement,
                                                                     make_plots=False,
                                                                     save_to_file=False,
                                                                     force_positive=force_positive,
                                                                     decay_mode="SPG",
                                                                     process_string="SPG",
                                                                     signs=["OS", "SS"],
                                                                     years=years,
                                                                     leptons=leptons,
                                                                     charges=[""],
                                                                     btags="",
                                                                     ptbins=ptbins)

        # replacement samples
        if options.replacement_samples:
            channelGeneratorReplacement.make_channel(lumi, extra_rebin=extra_rebin)
            channelGeneratorReplacement.make_channel(lumi, sign='OS', extra_rebin=extra_rebin)
            channelGeneratorReplacement.make_channel(lumi, sign='SS', extra_rebin=extra_rebin)

        # OS/SS plots
        if make_os_ss:
            for sign in signs:
                if not options.fit_only:
                    channelGenerator.make_channel(lumi, sign=sign, extra_rebin=extra_rebin, os_only=os_only)
                    if not options.inclusive_only:
                        for btag in btags:
                            channelGenerator.make_channel(lumi, sign=sign, btag=btag,
                                                          extra_rebin=extra_rebin * (btag_bin if btag != '0tag' else 1), os_only=os_only)
                if not options.inclusive_only:
                    for lepton in leptons:
                        if not options.fit_only:
                            channelGenerator.make_channel(lumi, sign=sign, lepton=lepton, extra_rebin=extra_rebin,
                                                          os_only=os_only)
                        for btag in btags:
                            if not options.fit_only:
                                channelGenerator.make_channel(lumi, sign=sign, btag=btag, lepton=lepton, extra_rebin=extra_rebin *
                                                              (btag_bin if btag != '0tag' else 1), os_only=os_only)
                            for charge in charges:
                                channelGenerator.make_channel(lumi, sign=sign, btag=btag, lepton=lepton, charge=charge,
                                                              extra_rebin=extra_rebin * (btag_bin if btag != '0tag' else 1), os_only=os_only)

        # OS-SS plots
        if not options.fit_only:
            if make_os_minus_ss:
                channelGenerator.make_channel(lumi, extra_rebin=extra_rebin, os_only=os_only)
                if not options.inclusive_only:
                    for year in years:
                        channelGenerator.make_channel([year], year=year, extra_rebin=extra_rebin, os_only=os_only)
                    for btag in btags:
                        channelGenerator.make_channel(lumi, btag=btag, extra_rebin=extra_rebin * (btag_bin if btag != '0tag' else 1),
                                                      os_only=os_only)
                    for lepton in leptons:
                        channelGenerator.make_channel(lumi, lepton=lepton, extra_rebin=extra_rebin, os_only=os_only)
                        for charge in charges:
                            channelGenerator.make_channel(lumi, lepton=lepton, charge=charge, extra_rebin=extra_rebin,
                                                          os_only=os_only)
                        for btag in btags:
                            channelGenerator.make_channel(lumi, btag=btag, lepton=lepton,
                                                          extra_rebin=extra_rebin * (btag_bin if btag != '0tag' else 1), os_only=os_only)
                            for charge in charges:
                                channelGenerator.make_channel(lumi, btag=btag, lepton=lepton, charge=charge,
                                                              extra_rebin=extra_rebin * (btag_bin if btag != '0tag' else 1), os_only=os_only)

        # add channels
        if options.replacement_samples:
            channelGenerator.get_config()['channels'].update(channelGeneratorReplacement.get_config()['channels'])

        out_name = f'{options.analysis_config}_{sample_config}{"_OSSS" if options.fit_type == "OS/SS" else ""}.yaml'
        with open(out_name, 'w') as outfile:
            yaml.dump(channelGenerator.get_config(), outfile, default_flow_style=False)


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-a', '--analysis-config',
                      action="store", dest="analysis_config",
                      help="analysis config name",
                      default="wplusd")
    parser.add_option('-d', '--decay-mode',
                      action="store", dest="decay_mode",
                      help="the decay mode string",
                      default="Dplus")
    parser.add_option('-s', '--samples',
                      action="store", dest="samples",
                      help="type of samples (truth, flavor, fit)",
                      default="truth")
    parser.add_option('-r', '--extra-rebin',
                      action="store", dest="extra_rebin",
                      help="extra rebin",
                      default=1)
    parser.add_option('-v', '--variables',
                      action="store", dest="variables",
                      help="the varaibels config (e.g. charmed_wjets, charmed_wjets_dstarPi0)",
                      default="charmed_wjets")
    parser.add_option('--samples-config',
                      action="store", dest="samples_config",
                      help="different sample configurations (madgraph_truth,sherpa_truth,sherpa2210_truth,madgraph_fxfx_truth)",
                      default="madgraph_truth")
    parser.add_option('--fit-only',
                      action="store_true", dest="fit_only",
                      help="only regions necessaty for the fit")
    parser.add_option('--fit-type',
                      action="store", dest="fit_type",
                      help="fit type (e.g. OS/SS or OS-SS)",
                      default="")
    parser.add_option('--inclusive-only',
                      action="store_true", dest="inclusive_only",
                      help="only inclusive regions (i.e. not split in charge / b-tag / etc..)")
    parser.add_option('--sys',
                      action="store_true", dest="systematics",
                      help="add systematics")
    parser.add_option('--replacement-samples',
                      action="store_true", dest="replacement_samples",
                      help="replace samples")
    parser.add_option('--process-string',
                      action="store", dest="process_string",
                      default="W#rightarrowl#nu+D, D#rightarrowK#pi#pi")

    # parse input arguments
    options, args = parser.parse_args()

    # run main
    main(options)
