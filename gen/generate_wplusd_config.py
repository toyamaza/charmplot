#!/usr/bin/env python
import gen.utils.templates as templates
import yaml

# process string
process_string = 'W#rightarrowl#nu+D, D#rightarrowK#pi#pi'


def main(options):

    for sample_config in options.samples.split(","):

        # TODO: make configurable?
        make_os_ss = True
        make_os_minus_ss = True
        force_positive = False

        # TODO: make configurable?
        samples = templates.WDTruthSamples().get()
        # samples = templates.WDFitSamples().get()

        # TODO: make configurable?
        signs = ['OS', 'SS']
        years = ['2017', '2018']
        leptons = ['el', 'mu']
        charges = ['plus', 'minus']
        btags = ['0tag', '1tag', '2tag']

        # systematics
        systematics = []
        if options.systematics:
            systematics = [
                'experimental',
                'matrix_method',
                'ttbar_theory_pdf',
                'ttbar_theory_choice',
                'ttbar_theory_qcd',
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
                channelGenerator.make_channel(years, sign=sign, extra_rebin=2)
                for btag in btags:
                    channelGenerator.make_channel(years, sign=sign, btag=btag, extra_rebin=2)
                for lepton in leptons:
                    channelGenerator.make_channel(years, sign=sign, lepton=lepton, extra_rebin=2)
                    for btag in btags:
                        channelGenerator.make_channel(years, sign=sign, btag=btag, lepton=lepton, extra_rebin=2)
                        for charge in charges:
                            channelGenerator.make_channel(years, sign=sign, btag=btag, lepton=lepton, charge=charge, extra_rebin=2)

        # OS-SS plots
        if make_os_minus_ss:
            channelGenerator.make_channel(years, extra_rebin=2)
            for btag in btags:
                channelGenerator.make_channel(years, btag=btag, extra_rebin=2)
            for lepton in leptons:
                channelGenerator.make_channel(years, lepton=lepton, extra_rebin=2)
                for btag in btags:
                    channelGenerator.make_channel(years, btag=btag, lepton=lepton, extra_rebin=2)
                    for charge in charges:
                        channelGenerator.make_channel(years, btag=btag, lepton=lepton, charge=charge, extra_rebin=2)

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
                      default="SR_Dplus")
    parser.add_option('-s', '--samples',
                      action="store", dest="samples",
                      help="different sample configurations (madgraph_truth,sherpa_truth,sherpa2210_truth,madgraph_fxfx_truth)",
                      default="madgraph_truth")
    parser.add_option('--sys',
                      action="store_true", dest="systematics",
                      help="add systematics")

    # parse input arguments
    options, args = parser.parse_args()

    # run main
    main(options)
