#!/usr/bin/env python
import gen.utils.templates as templates
import sys
import yaml


def main(options):

    # additional rebinning
    extra_rebin = float(options.extra_rebin)

    # TODO: make configurable?
    make_os_ss = True
    make_os_minus_ss = not options.fit_only
    os_only = options.fit_only
    force_positive = True
    if options.fit_type == "OS-SS":
        force_positive = False

    # sample type
    if options.samples.lower() == 'truth':
        samples = templates.WDTruthSamples(fitType=options.fit_type,
                                           MockMC=options.fit_type != '',
                                           decayMode=options.decay_mode,
                                           truthDiffBins=options.truth_differential_bins,
                                           splitSignalSamples=options.split_signal_samples)
        if options.replacement_samples:
            samples_replacement = templates.ReplacementSamples(truthDiffBins=options.truth_differential_bins,
                                                               splitSignalSamples=options.split_signal_samples)
    elif options.samples.lower() == 'flavor' or options.samples.lower() == 'flavour':
        samples = templates.WDFlavourSamples()
    elif options.samples.lower() == 'fit':
        samples = templates.WDFitSamples()
    elif options.samples.lower() == 'spg_comparison':
        samples = templates.SPGComparison(truthDiffBins=options.truth_differential_bins,
                                          splitSignalSamples=options.split_signal_samples, decay_mode=options.decay_mode)
    elif options.samples.lower() == 'bkg_comparison':
        samples = templates.BKGComparison()
    elif options.samples.lower() == 'wplusd_comparison':
        samples = templates.WDComparisonSamples()
    else:
        print(f"ERROR: Unknown samples type {options.samples}")
        sys.exit(1)

    # TODO: make configurable?
    signs = ['OS', 'SS']
    lumi = ['2015', '2016A-B', '2016C-L', '2017', '2018']
    years = []  # to be used without 'split_by_period' option in charmpp
    leptons = ['mu', 'el']
    charges = ['minus', 'plus']
    btags = ['0tag', '1tag']
    ptbins = ['']
    if options.differential_bins:
        ptbins = ['pt_bin1', 'pt_bin2', 'pt_bin3', 'pt_bin4', 'pt_bin5', '']

    # binning for the 1-tag region
    # use '40' for 1 bin for D+
    btag_bin = 40

    # replace samples to increase mc stats
    replacement_samples = {}
    if options.replacement_samples:
        if options.truth_differential_bins:
            replacement_samples = {
                'Wjets_emu_Matched_truth_pt_bin1': '<charge>_Replacement_Matched_truth_pt_bin1',
                'Wjets_emu_Matched_truth_pt_bin2': '<charge>_Replacement_Matched_truth_pt_bin2',
                'Wjets_emu_Matched_truth_pt_bin3': '<charge>_Replacement_Matched_truth_pt_bin3',
                'Wjets_emu_Matched_truth_pt_bin4': '<charge>_Replacement_Matched_truth_pt_bin4',
                'Wjets_emu_Matched_truth_pt_bin5': '<charge>_Replacement_Matched_truth_pt_bin5',
            }
        else:
            replacement_samples = {'Wjets_emu_Matched': '<charge>_Replacement_Matched'}
        replacement_samples.update({
            'Wjets_emu_411MisMatched': '<charge>_Replacement_411MisMatched',
            'Wjets_emu_Charm': '<charge>_Replacement_CharmMisMatched',
            'Wjets_emu_Rest': '<charge>_Replacement_Wjets_emu_Rest',
            'Wjets_emu_MisMatched': '<charge>_Replacement_Wjets_emu_MisMatched',
            'DibosonVjetsTau': '<charge>_Replacement_DibosonVjetsTau',
        })

    # systematics
    systematics = []
    if options.systematics:
        systematics = [
            'experimental',
            'matrix_method',
            'ttbar_theory_alt_samples',
            'ttbar_theory_choice',
            'ttbar_theory_pdf',
            'ttbar_theory_qcd',
            'wjets_theory',
            # 'proxy_norm',
        ]

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
                                                  ptbins=ptbins,
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
    parser.add_option('--differential-bins',
                      action="store_true", dest="differential_bins",
                      default=False)
    parser.add_option('--truth-differential-bins',
                      action="store_true", dest="truth_differential_bins",
                      default=False)
    parser.add_option('--split-signal-samples',
                      action="store_true", dest="split_signal_samples",
                      default=False)

    # parse input arguments
    options, args = parser.parse_args()

    # run main
    main(options)
