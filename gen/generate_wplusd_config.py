#!/usr/bin/env python
import gen.utils.templates as templates
import sys
import yaml


def main(options):

    # additional rebinning
    extra_rebin = float(options.extra_rebin)

    # TODO: make configurable?
    make_os_ss = not options.skip_os
    make_os_minus_ss = not options.fit_only and not options.skip_os_ss
    force_positive = True
    if options.not_force_positive:
        force_positive = False

    # sample type
    if options.samples.lower() == 'truth':
        samples = templates.WDTruthSamples(fitType=options.fit_type,
                                           MockMC=options.fit_type != '',
                                           decayMode=options.decay_mode,
                                           truthDiffBins=options.truth_differential_bins,
                                           samplesConfOverride=options.samples_config,
                                           eta_bins=options.differential_eta)
    elif options.samples.lower() == 'flavor' or options.samples.lower() == 'flavour':
        samples = templates.WDFlavourSamples()
    elif options.samples.lower() == 'fit':
        samples = templates.WDFitSamples()
    elif options.samples.lower() == 'fake_track' or options.samples.lower() == 'fake_track':
        samples = templates.WDFakeTrackSamples()
    elif options.samples.lower() == 'spg_comparison':
        samples = templates.SPGComparison(truthDiffBins=options.truth_differential_bins,
                                          decay_mode=options.decay_mode,
                                          signal_only=options.spg_signal_only,
                                          background_only=options.spg_background_only,
                                          eta_bins=options.differential_eta)
    elif options.samples.lower() == 'bkg_comparison':
        samples = templates.BKGComparison(decay_mode=options.decay_mode)
    elif options.samples.lower() == 'signal_comparison':
        samples = templates.SignalComparison(decay_mode=options.decay_mode)
    elif options.samples.lower() == 'wplusd_comparison':
        samples = templates.WDComparisonSamples()
    elif options.samples.lower() == 'wjets_sherpa_sys':
        samples = templates.WjetsSherpaSys(decay_mode=options.decay_mode)
    elif options.samples.lower() == 'spg_sys':
        samples = templates.SPGSysComparison()
    else:
        print(f"ERROR: Unknown samples type {options.samples}")
        sys.exit(1)

    # replacement samples
    if options.replacement_samples:
        samples_replacement = templates.ReplacementSamples(truthDiffBins=options.truth_differential_bins,
                                                           decay_mode=options.decay_mode,
                                                           comparison=options.samples.lower(),
                                                           eta_bins=options.differential_eta)

    # TODO: make configurable?
    signs = ['OS', 'SS']
    # lumi = ['2015', '2016A-B', '2016C-L', '2017', '2018']
    lumi = ['2015', '2016', '2017', '2018']
    years = []  # to be used without 'split_by_period' option in charmpp
    leptons = ['mu', 'el']
    charges = ['minus', 'plus']
    btags = ['0tag', '1tag']
    differential_bins = ['']
    if options.differential_bins:
        if options.differential_eta:
            differential_bins = ['eta_bin1', 'eta_bin2', 'eta_bin3', 'eta_bin4', 'eta_bin5', '']
        else:
            differential_bins = ['pt_bin1', 'pt_bin2', 'pt_bin3', 'pt_bin4', 'pt_bin5', '']

    # override from CLI
    if options.btags:
        btags = options.btags.split(",")
        if type(btags) != list:
            btags = [btags]

    if options.leptons:
        leptons = options.leptons.split(",")
        if type(leptons) != list:
            leptons = [leptons]

    if options.charges:
        charges = options.charges.split(",")
        if type(charges) != list:
            charges = [charges]

    # binning for the 1-tag region
    # use '-1' for 1 bin for D+
    btag_bin = -1

    # replace samples to increase mc stats
    replacement_samples = {}
    if options.replacement_samples:
        if options.truth_differential_bins:
            if options.decay_mode == "Dplus":
                pass
            else:
                if options.differential_eta:
                    replacement_samples = {
                        'Sherpa2211_WplusD_Matched_truth_eta_bin1': '<sign>_Replacement_MatchedMG<charge>_truth_eta_bin1',
                        'Sherpa2211_WplusD_Matched_truth_eta_bin2': '<sign>_Replacement_MatchedMG<charge>_truth_eta_bin2',
                        'Sherpa2211_WplusD_Matched_truth_eta_bin3': '<sign>_Replacement_MatchedMG<charge>_truth_eta_bin3',
                        'Sherpa2211_WplusD_Matched_truth_eta_bin4': '<sign>_Replacement_MatchedMG<charge>_truth_eta_bin4',
                        'Sherpa2211_WplusD_Matched_truth_eta_bin5': '<sign>_Replacement_MatchedMG<charge>_truth_eta_bin5',
                    }
                else:
                    replacement_samples = {
                        'Sherpa2211_WplusD_Matched_truth_pt_bin1': '<sign>_Replacement_MatchedMG<charge>_truth_pt_bin1',
                        'Sherpa2211_WplusD_Matched_truth_pt_bin2': '<sign>_Replacement_MatchedMG<charge>_truth_pt_bin2',
                        'Sherpa2211_WplusD_Matched_truth_pt_bin3': '<sign>_Replacement_MatchedMG<charge>_truth_pt_bin3',
                        'Sherpa2211_WplusD_Matched_truth_pt_bin4': '<sign>_Replacement_MatchedMG<charge>_truth_pt_bin4',
                        'Sherpa2211_WplusD_Matched_truth_pt_bin5': '<sign>_Replacement_MatchedMG<charge>_truth_pt_bin5',
                    }
        else:
            if options.decay_mode == "Dplus":
                pass
            else:
                replacement_samples = {
                    'Sherpa2211_WplusD_Matched': '<sign>_Replacement_MatchedMG<charge>',
                }
        replacement_samples.update({
            # 'MG_Wjets_Charm': '<sign>_Replacement_CharmMisMatched',
            'MG_Wjets_MisMatched': '<sign>_Replacement_MisMatched',
            'Sherpa2211_Wjets_MisMatched': '<sign>_Replacement_MisMatched',
            'MG_Wjets_Rest': '<sign>_Replacement_Rest',
            'Sherpa2211_Wjets_Rest': '<sign>_Replacement_Rest',
        })

    # systematics
    systematics = []
    if options.systematics:
        systematics = [
            'experimental',
            'matrix_method',
            'theory_prod_frac',
            'ttbar_theory_alt_samples',
            'ttbar_theory_choice',
            'ttbar_theory_pdf',
            'ttbar_theory_qcd',
            'wjets_theory',
            'wjets_theory_madgraph',
            'sherpa2211_theory_qcd_fit',
            'sherpa2211_theory_as_fit',
            'sherpa2211_theory_pdf',
            'wjets_mismatch_alt_samples',
            'dplus_resolution_signal',
            'dplus_resolution_background',
            'dplus_signal_morphing',
        ]
        if options.decay_mode == "Dplus":
            systematics += [
                'sherpa2211_wjets_bkg_alt_samples',
            ]
            if options.truth_differential_bins:
                if options.differential_eta:
                    systematics += [
                        'wplusd_signal_alt_samples_eta_bin1',
                        'wplusd_signal_alt_samples_eta_bin2',
                        'wplusd_signal_alt_samples_eta_bin3',
                        'wplusd_signal_alt_samples_eta_bin4',
                        'wplusd_signal_alt_samples_eta_bin5',
                    ]
                else:
                    systematics += [
                        'wplusd_signal_alt_samples_pt_bin1',
                        'wplusd_signal_alt_samples_pt_bin2',
                        'wplusd_signal_alt_samples_pt_bin3',
                        'wplusd_signal_alt_samples_pt_bin4',
                        'wplusd_signal_alt_samples_pt_bin5',
                    ]
            else:
                systematics += [
                    'wplusd_signal_alt_samples',
                ]
        elif options.decay_mode == "Dstar":
            systematics += [
                'wjets_bkg_alt_samples',
                'wjets_bkg_alt_samples_1tag',
                'wjets_bkg_fit',
            ]

    if options.sys_configs:
        systematics = options.sys_configs.split(",")

    # Base config
    config = templates.DataMCConfig(variables=options.variables,
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
                                                  differential_bins=differential_bins,
                                                  process_string=options.process_string,
                                                  force_positive=force_positive,
                                                  replacement_samples=replacement_samples,
                                                  os_ss_sub=(options.fit_type == "OS-SS"))

    # Helper object to generate channels
    if options.replacement_samples:
        channelGeneratorReplacement = templates.ChannelGenerator(config=config,
                                                                 samples=samples_replacement,
                                                                 make_plots=False,
                                                                 save_to_file=False,
                                                                 force_positive=force_positive,
                                                                 decay_mode="Replacement",
                                                                 process_string="Replacement",
                                                                 signs=["OS", "SS"],
                                                                 years=years,
                                                                 leptons=leptons,
                                                                 charges=[""],
                                                                 btags="",
                                                                 differential_bins=differential_bins)

    # replacement samples
    if options.replacement_samples:
        channelGeneratorReplacement.make_channel(lumi, extra_rebin=extra_rebin)
        channelGeneratorReplacement.make_channel(lumi, sign='OS', extra_rebin=extra_rebin)
        channelGeneratorReplacement.make_channel(lumi, sign='SS', extra_rebin=extra_rebin)

    # OS/SS plots
    if make_os_ss:
        for sign in signs:
            if not options.fit_only:
                if len(btags) > 1:
                    channelGenerator.make_channel(lumi, sign=sign, extra_rebin=extra_rebin)
                for btag in btags:
                    channelGenerator.make_channel(lumi, sign=sign, btag=btag,
                                                  extra_rebin=extra_rebin * (btag_bin if btag != '0tag' else 1))
            if not options.inclusive_only:
                for lepton in leptons:
                    if not options.fit_only:
                        channelGenerator.make_channel(lumi, sign=sign, lepton=lepton, extra_rebin=extra_rebin)
                    for btag in btags:
                        if not options.fit_only:
                            channelGenerator.make_channel(lumi, sign=sign, btag=btag, lepton=lepton, extra_rebin=extra_rebin *
                                                          (btag_bin if btag != '0tag' else 1))
                        for charge in charges:
                            channelGenerator.make_channel(lumi, sign=sign, btag=btag, lepton=lepton, charge=charge,
                                                          extra_rebin=extra_rebin * (btag_bin if btag != '0tag' else 1))

            if not options.skip_inclusive_leptons:
                for charge in charges:
                    if not options.fit_only:
                        channelGenerator.make_channel(lumi, sign=sign, charge=charge,
                                                      extra_rebin=extra_rebin)
                    for btag in btags:
                        channelGenerator.make_channel(lumi, sign=sign, btag=btag, charge=charge,
                                                      extra_rebin=extra_rebin * (btag_bin if btag != '0tag' else 1))

    # OS-SS plots
    if not options.fit_only:
        if make_os_minus_ss:
            if len(btags) > 1:
                channelGenerator.make_channel(lumi, extra_rebin=extra_rebin)
                for btag in btags:
                    channelGenerator.make_channel(lumi, btag=btag,
                                                  extra_rebin=extra_rebin * (btag_bin if btag != '0tag' else 1))
            else:
                for btag in btags:
                    channelGenerator.make_channel(lumi, btag=btag, extra_rebin=extra_rebin * (btag_bin if btag != '0tag' else 1))
            if not options.inclusive_only:
                for year in years:
                    channelGenerator.make_channel([year], year=year, extra_rebin=extra_rebin)
                for btag in btags:
                    channelGenerator.make_channel(lumi, btag=btag, extra_rebin=extra_rebin * (btag_bin if btag != '0tag' else 1))
                for lepton in leptons:
                    channelGenerator.make_channel(lumi, lepton=lepton, extra_rebin=extra_rebin)
                    for charge in charges:
                        channelGenerator.make_channel(lumi, lepton=lepton, charge=charge, extra_rebin=extra_rebin)
                    for btag in btags:
                        channelGenerator.make_channel(lumi, btag=btag, lepton=lepton,
                                                      extra_rebin=extra_rebin * (btag_bin if btag != '0tag' else 1))
                        for charge in charges:
                            channelGenerator.make_channel(lumi, btag=btag, lepton=lepton, charge=charge,
                                                          extra_rebin=extra_rebin * (btag_bin if btag != '0tag' else 1))

            if not options.skip_inclusive_leptons:
                for charge in charges:
                    channelGenerator.make_channel(lumi, charge=charge,
                                                  extra_rebin=extra_rebin)
                    for btag in btags:
                        channelGenerator.make_channel(lumi, btag=btag, charge=charge,
                                                      extra_rebin=extra_rebin * (btag_bin if btag != '0tag' else 1))

    # add channels
    if options.replacement_samples:
        channelGenerator.get_config()['channels'].update(channelGeneratorReplacement.get_config()['channels'])

    out_name = f'{options.analysis_config}{"_OSSS" if options.fit_type == "OS/SS" else ""}.yaml'
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
                      help="type of samples (truth, flavor, fit, fake_track)",
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
                      default="")
    parser.add_option('--fit-only',
                      action="store_true", dest="fit_only",
                      help="only regions necessaty for the fit")
    parser.add_option('--skip-os-ss',
                      action="store_true", dest="skip_os_ss",
                      default=False, help="don't make OS-SS plots")
    parser.add_option('--skip-os',
                      action="store_true", dest="skip_os",
                      default=False, help="don't make OS / SS plots")
    parser.add_option('--skip-inclusive-leptons',
                      action="store_true", dest="skip_inclusive_leptons",
                      default=False, help="don't make inclusive W+ / W- plots")
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
    parser.add_option('--sys-configs',
                      action="store", dest="sys_configs",
                      help="comma sepparated list of sys configs")
    parser.add_option('--replacement-samples',
                      action="store_true", dest="replacement_samples",
                      help="replace samples")
    parser.add_option('--process-string',
                      action="store", dest="process_string",
                      default="W#rightarrowl#nu+D, D#rightarrowK#pi#pi")
    parser.add_option('--differential-bins',
                      action="store_true", dest="differential_bins",
                      default=False)
    parser.add_option('--truth-differential-bins',
                      action="store_true", dest="truth_differential_bins",
                      default=False)
    parser.add_option('--differential-eta',
                      action="store_true", dest="differential_eta",
                      default=False)
    parser.add_option('--spg-signal-only',
                      action="store_true", dest="spg_signal_only",
                      default=False)
    parser.add_option('--spg-background-only',
                      action="store_true", dest="spg_background_only",
                      default=False)
    parser.add_option('--not-force-positive',
                      action="store_true", dest="not_force_positive",
                      default=False)
    # ----------------------------------------------------
    # arguments: regions / channels override
    # ----------------------------------------------------
    parser.add_option('--btags',
                      action="store", dest="btags",
                      default="")
    parser.add_option('--leptons',
                      action="store", dest="leptons",
                      default="")
    parser.add_option('--charges',
                      action="store", dest="charges",
                      default="")

    # parse input arguments
    options, args = parser.parse_args()

    # run main
    main(options)
