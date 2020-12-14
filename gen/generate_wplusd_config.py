#!/usr/bin/env python
import gen.utils.templates as templates
import sys
import yaml

# process string
process_string = 'W#rightarrowl#nu+D, D#rightarrowK#pi#pi'


def main(options):

    for sample_config in options.samples_config.split(","):

        # TODO: make configurable
        extra_rebin = float(options.extra_rebin)

        # TODO: make configurable?
        make_os_ss = True
        make_os_minus_ss = not options.fit_only
        force_positive = options.samples.lower() == 'fit'

        # sample type
        if options.samples.lower() == 'truth':
            samples = templates.WDTruthSamples().get()
        elif options.samples.lower() == 'flavor' or options.samples.lower() == 'flavour':
            samples = templates.WDFlavourSamples().get()
        elif options.samples.lower() == 'fit':
            samples = templates.WDFitSamples().get()
        else:
            print(f"ERROR: Unknown samples type {options.samples}")
            sys.exit(1)

        # TODO: make configurable?
        signs = ['OS', 'SS']
        years = ['2017', '2018']
        leptons = ['el', 'mu']
        charges = ['plus', 'minus']
        btags = ['0tag', ['1tag', '2tag']]

        # systematics
        systematics = []
        if options.systematics:
            systematics = [
                'experimental',
                'matrix_method',
                'ttbar_theory_pdf',
                'ttbar_theory_choice',
                'ttbar_theory_qcd',
                'ttbar_theory_alt_samples',
                'wjets_theory_alt_samples',
            ]

        # Base config
        config = templates.DataMCConfig(variables='charmed_wjets',
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
                                                      process_string=process_string,
                                                      sample_config=sample_config,
                                                      force_positive=force_positive)

        # OS/SS plots
        if make_os_ss:
            for sign in signs:
                if not options.fit_only:
                    channelGenerator.make_channel(years, sign=sign, extra_rebin=extra_rebin)
                    for btag in btags:
                        channelGenerator.make_channel(years, sign=sign, btag=btag, extra_rebin=extra_rebin * (2 if btag != '0tag' else 1))
                for lepton in leptons:
                    if not options.fit_only:
                        channelGenerator.make_channel(years, sign=sign, lepton=lepton, extra_rebin=extra_rebin)
                    for btag in btags:
                        if not options.fit_only:
                            channelGenerator.make_channel(years, sign=sign, btag=btag, lepton=lepton, extra_rebin=extra_rebin * (2 if btag != '0tag' else 1))
                        for charge in charges:
                            channelGenerator.make_channel(years, sign=sign, btag=btag, lepton=lepton, charge=charge, extra_rebin=extra_rebin * (2 if btag != '0tag' else 1))

        if not options.fit_only:
            # OS-SS plots
            if make_os_minus_ss:
                channelGenerator.make_channel(years, extra_rebin=extra_rebin)
                for btag in btags:
                    channelGenerator.make_channel(years, btag=btag, extra_rebin=extra_rebin * (2 if btag != '0tag' else 1))
                for lepton in leptons:
                    channelGenerator.make_channel(years, lepton=lepton, extra_rebin=extra_rebin)
                    for charge in charges:
                        channelGenerator.make_channel(years, lepton=lepton, charge=charge, extra_rebin=extra_rebin)
                    for btag in btags:
                        channelGenerator.make_channel(years, btag=btag, lepton=lepton, extra_rebin=extra_rebin * (2 if btag != '0tag' else 1))
                        for charge in charges:
                            channelGenerator.make_channel(years, btag=btag, lepton=lepton, charge=charge, extra_rebin=extra_rebin * (2 if btag != '0tag' else 1))

        with open(f'{options.analysis_config}_{sample_config}.yaml', 'w') as outfile:
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
                      default="wplud")
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
    parser.add_option('--samples-config',
                      action="store", dest="samples_config",
                      help="different sample configurations (madgraph_truth,sherpa_truth,sherpa2210_truth,madgraph_fxfx_truth)",
                      default="madgraph_truth")
    parser.add_option('--fit-only',
                      action="store_true", dest="fit_only",
                      help="only regions necessaty for the fit")
    parser.add_option('--sys',
                      action="store_true", dest="systematics",
                      help="add systematics")

    # parse input arguments
    options, args = parser.parse_args()

    # run main
    main(options)
