#!/usr/bin/env python
import gen.utils.proxies as proxies


def flatten(xs):
    result = []
    if isinstance(xs, (list, tuple)):
        for x in xs:
            result.extend(flatten(x))
    else:
        result.append(xs)
    return result


label_dict = {
    'madgraph_truth': 'MadGraph LO',
    'sherpa_truth': 'Sherpa 2.2.1',
    'sherpa2210_truth': 'Sherpa 2.2.10',
    'madgraph_fxfx_truth': 'MadGraph FxFx',
    'madgraph': 'MadGraph LO',
    'sherpa': 'Sherpa 2.2.1',
    'sherpa2210': 'Sherpa 2.2.10',
    'madgraph_fxfx': 'MadGraph FxFx',
    'truth_comparison': 'LO MG vs. LO Powheg',
    'wplusd_comparison': 'MadGraph LO',
    'flavor_comparison': 'MadGraph LO',
    'spg_comparison': 'MadGraph LO vs. SPG',
    'multijet_comparison': 'Matrix Method vs Fake Factor',
    'multijet_composition': 'Matrix Method',
}


def format(string):
    if string:
        return f'{string}_'
    else:
        return string


class DataMCConfig:

    def __init__(self, variables, sample_config, systematics=[]):
        self.variables = variables
        self.sample_config = sample_config
        self.systematics = systematics

    def to_dict(self):
        out = {
            'variablesConf': self.variables,
            'samplesConf': self.sample_config,
            'channels': {},
        }
        if self.sample_config not in ['truth_comparison', 'wplusd_comparison', 'spg_comparison', 'multijet_comparison', 'multijet_composition']:
            out['data'] = 'Data'
        if self.systematics:
            out['systematics'] = self.systematics
        return out


class WDFitSamples:

    def __init__(self, loose_sr=False):
        self.loose_sr = loose_sr
        self.samples = [
            ['MockMC', proxies.MockMC()],
            ['Wjets_emu_Matched', proxies.Matched(os_minus_ss_fit_configuration=True)],
            ['Wjets_cjets_emu_Rest', proxies.Rest(os_minus_ss_fit_configuration=True, loose_sr=self.loose_sr)],
            ['Wjets_bujets_emu_Rest', proxies.Rest(os_minus_ss_fit_configuration=True, loose_sr=self.loose_sr)],
            ['Top_Matched', proxies.Matched(os_minus_ss_fit_configuration=True)],
            ['Top_Rest', proxies.Rest(os_minus_ss_fit_configuration=True)],
            ['Wjets_emu_NoMatch', proxies.NoMatch(os_minus_ss_fit_configuration=True)],
            ['Other', proxies.PlainChannel(os_minus_ss_fit_configuration=True)],
            ['Zjets_emu', proxies.PlainChannel(os_minus_ss_fit_configuration=True)],
            ['Multijet_MatrixMethod', proxies.MatrixMethod(os_minus_ss_fit_configuration=True)]
        ]

    def get(self):
        return self.samples


class WDTruthSamplesNew:

    def __init__(self, os_minus_ss_fit_configuration=False, loose_sr=False, OS_and_SS_fit=False, MockMC=True, decayMode="Dplus"):
        self.os_minus_ss_fit_configuration = os_minus_ss_fit_configuration
        self.loose_sr = loose_sr
        self.OS_and_SS_fit = OS_and_SS_fit
        self.MockMC = MockMC
        self.decayMode = decayMode
        self.samples = []
        if self.os_minus_ss_fit_configuration and self.MockMC:
            self.samples += [['MockMC', proxies.MockMC(subtract_mj=False)]]
        self.samples += [
            ['Wjets_emu_Matched', proxies.Matched(os_minus_ss_fit_configuration=self.os_minus_ss_fit_configuration)],
        ]
        if self.decayMode == "Dplus":
            self.samples += [
                ['Wjets_emu_411MisMatched', proxies.MisMatched(
                    os_minus_ss_fit_configuration=self.os_minus_ss_fit_configuration, loose_sr=self.loose_sr, pdgId="411")],
            ]
        elif self.decayMode == "DstarKPiPi0":
            self.samples += [
                ['Wjets_emu_413MisMatched', proxies.MisMatched(
                    os_minus_ss_fit_configuration=self.os_minus_ss_fit_configuration, loose_sr=self.loose_sr, pdgId="413")],
            ]
        self.samples += [
            ['Wjets_emu_Charm', proxies.MatchedCharm(os_minus_ss_fit_configuration=self.os_minus_ss_fit_configuration, loose_sr=self.loose_sr, decayMode=self.decayMode)],  # noqa: E501
            ['Wjets_emu_MisMatched', proxies.GenericChannel(name="MisMatched", os_minus_ss_fit_configuration=self.os_minus_ss_fit_configuration, region=["MisMatched", "MatchedNoFid"])],  # noqa: E501
            ['Wjets_Rest', proxies.NoMatchBackground(os_minus_ss_fit_configuration=self.os_minus_ss_fit_configuration, add_no_truth_match=True, loose_sr=self.loose_sr)],  # noqa: E501
            ['Top', proxies.PlainChannel(os_minus_ss_fit_configuration=self.os_minus_ss_fit_configuration)],
            ['DibosonZjetsTau', proxies.PlainChannel(os_minus_ss_fit_configuration=self.os_minus_ss_fit_configuration)],
        ]
        if self.os_minus_ss_fit_configuration:
            self.samples += [['Multijet_MatrixMethod', proxies.MatrixMethod(os_minus_ss_fit_configuration=True, fake_factor=False)]]
        else:
            self.samples += [['Multijet_MatrixMethod', proxies.MatrixMethod(fake_factor=False)]]
            if self.MockMC:
                self.samples += [['MockMC_minus_MC', proxies.MockMC(subtract_mj=True)]]

    def get(self):
        return self.samples


class MultiJetComparison:

    samples = [
        ['Multijet_MatrixMethod', proxies.MatrixMethod(fake_factor=False)],
        ['Multijet_MatrixMethod_FF', proxies.MatrixMethod(fake_factor=True)],
    ]

    def get(self):
        return self.samples


class MultiJetComposition:

    samples = [
        ['Multijet_MatrixMethod_AntiTight', proxies.MatrixMethod(loose_only=True)],
        ['Multijet_MatrixMethod_Tight', proxies.MatrixMethod(tight_only=True)],
        ['Multijet_MatrixMethod', proxies.MatrixMethod()],
    ]

    def get(self):
        return self.samples


class WDTruthSamples:

    samples = [
        ['Wjets_emu_Matched', proxies.Matched()],
        ['Wjets_cjets_emu_Rest', proxies.Rest()],
        ['Wjets_bujets_emu_Rest', proxies.Rest()],
        ['Top_Matched', proxies.Matched()],
        ['Top_Rest', proxies.Rest()],
        ['Wjets_emu_NoMatch', proxies.NoMatch()],
        ['Other'],
        ['Zjets_emu'],
        ['Multijet_MatrixMethod', proxies.MatrixMethod()],
    ]

    def get(self):
        return self.samples


class WDFlavourSamples:

    samples = [
        ['Wjets_cjets_emu'],
        ['Wjets_bjets_emu'],
        ['Wjets_light_emu'],
        ['Top'],
        ['Zjets_emu'],
        ['Other'],
        ['Multijet_MatrixMethod', proxies.MatrixMethod()]
    ]

    def get(self):
        return self.samples


class WDTruthComparisonSamples:

    samples = {
        'Matched': [
            ['Wjets_emu_Matched_OS-SS', proxies.GenericChannel(region="Matched", name="MatchedOS-SS", os_minus_ss_fit_configuration=True)],
            # ['Wjets_emu_Matched_OS', proxies.GenericChannel(region="Matched", name="MatchedOS")],
            # ['Wjets_emu_Matched_SS', proxies.GenericChannel(region="Matched", name="MatchedSS", os_minus_ss_fit_configuration=True, ss_only=True)],
            ['SPG_Matched', proxies.GenericChannel(name="MatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS_Matched", "-inclusive_Dplus_SS_Matched"])],
        ],
        'Matched_OppositeSign': [
            ['Wjets_emu_Matched_OS', proxies.GenericChannel(region="Matched", name="MatchedOS")],
            ['SPG_Matched', proxies.GenericChannel(name="MatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS_Matched", "-inclusive_Dplus_SS_Matched"])],
        ],
        'Matched_SameSign': [
            ['Wjets_emu_Matched_SS', proxies.GenericChannel(region="Matched", name="MatchedSS", os_minus_ss_fit_configuration=True, ss_only=True)],
            ['SPG_Matched', proxies.GenericChannel(name="MatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS_Matched", "-inclusive_Dplus_SS_Matched"])],
        ],
        '411MisMatched': [
            ['Wjets_emu_411MisMatched_OS-SS', proxies.GenericChannel(region="411MisMatched", name="411MisMatchedOS-SS", os_minus_ss_fit_configuration=True)],
            # ['Wjets_emu_411MisMatched_OS', proxies.GenericChannel(region="411MisMatched", name="411MisMatchedOS")],
            # ['Wjets_emu_411MisMatched_SS', proxies.GenericChannel(region="411MisMatched", name="411MisMatchedSS", os_minus_ss_fit_configuration=True, ss_only=True)],
            ['SPG_411MisMatched', proxies.GenericChannel(name="411MisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS_411MisMatched", "-inclusive_Dplus_SS_411MisMatched"])],  # noqa: E501
        ],
        '411MisMatched_OppositeSign': [
            ['Wjets_emu_411MisMatched_OS', proxies.GenericChannel(region="411MisMatched", name="411MisMatchedOS")],
            ['SPG_411MisMatched', proxies.GenericChannel(name="411MisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS_411MisMatched", "-inclusive_Dplus_SS_411MisMatched"])],  # noqa: E501
        ],
        '411MisMatched_SameSign': [
            ['Wjets_emu_411MisMatched_SS', proxies.GenericChannel(region="411MisMatched", name="411MisMatchedSS", os_minus_ss_fit_configuration=True, ss_only=True)],
            ['SPG_411MisMatched', proxies.GenericChannel(name="411MisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS_411MisMatched", "-inclusive_Dplus_SS_411MisMatched"])],  # noqa: E501
        ],
        '431MisMatched': [
            ['Wjets_emu_431MisMatched_OS-SS', proxies.GenericChannel(region="431MisMatched", name="431MisMatchedOS-SS", os_minus_ss_fit_configuration=True)],
            # ['Wjets_emu_431MisMatched_OS', proxies.GenericChannel(region="431MisMatched", name="431MisMatchedOS")],
            # ['Wjets_emu_431MisMatched_SS', proxies.GenericChannel(region="431MisMatched", name="431MisMatchedSS", os_minus_ss_fit_configuration=True, ss_only=True)],
            ['SPG_431MisMatched', proxies.GenericChannel(name="431MisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS", "-inclusive_Dplus_SS"])],
        ],
        '431MisMatched_OppositeSign': [
            ['Wjets_emu_431MisMatched_OS', proxies.GenericChannel(region="431MisMatched", name="431MisMatchedOS")],
            ['SPG_431MisMatched', proxies.GenericChannel(name="431MisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS"])],
            # ['SPG_431MisMatched', proxies.GenericChannel(name="431MisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS", "-inclusive_Dplus_SS"])],
        ],
        '431MisMatched_SameSign': [
            ['Wjets_emu_431MisMatched_SS', proxies.GenericChannel(region="431MisMatched", name="431MisMatchedSS", os_minus_ss_fit_configuration=True, ss_only=True)],
            ['SPG_431MisMatched', proxies.GenericChannel(name="431MisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_SS"])],
            # ['SPG_431MisMatched', proxies.GenericChannel(name="431MisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS", "-inclusive_Dplus_SS"])],
        ],
        '421MisMatched': [
            ['Wjets_emu_421MisMatched_OS-SS', proxies.GenericChannel(region=["421MisMatched", "413MisMatched"], name="421MisMatchedOS-SS", os_minus_ss_fit_configuration=True)],  # noqa: E501
            # ['Wjets_emu_421MisMatched_OS', proxies.GenericChannel(region=["421MisMatched", "413MisMatched"], name="421MisMatchedOS")],
            # ['Wjets_emu_421MisMatched_SS', proxies.GenericChannel(region=["421MisMatched", "413MisMatched"], name="421MisMatchedSS", os_minus_ss_fit_configuration=True, ss_only=True)],
            ['SPG_421MisMatched', proxies.GenericChannel(name="421MisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS", "-inclusive_Dplus_SS"])],
        ],
        '421MisMatched_OppositeSign': [
            ['Wjets_emu_421MisMatched_OS', proxies.GenericChannel(region=["421MisMatched", "413MisMatched"], name="421MisMatchedOS")],
            ['SPG_421MisMatched', proxies.GenericChannel(name="421MisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS"])],
        ],
        '421MisMatched_SameSign': [
            ['Wjets_emu_421MisMatched_SS', proxies.GenericChannel(region=["421MisMatched", "413MisMatched"], name="421MisMatchedSS", os_minus_ss_fit_configuration=True, ss_only=True)],
            ['SPG_421MisMatched', proxies.GenericChannel(name="421MisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_SS"])],
        ],
        'BaryonMisMatched': [
            ['Wjets_emu_BaryonMisMatched_OS-SS', proxies.GenericChannel(region="BaryonMisMatched", name="BaryonMisMatchedOS-SS", os_minus_ss_fit_configuration=True)],  # noqa: E501
            # ['Wjets_emu_BaryonMisMatched_OS', proxies.GenericChannel(region="BaryonMisMatched", name="BaryonMisMatchedOS")],
            # ['Wjets_emu_BaryonMisMatched_SS', proxies.GenericChannel(region="BaryonMisMatched", name="BaryonMisMatchedSS", os_minus_ss_fit_configuration=True, ss_only=True)],
            ['SPG_BaryonMisMatched', proxies.GenericChannel(name="BaryonMisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS", "-inclusive_Dplus_SS"])],  # noqa: E501
        ],
        'BaryonMisMatched_OppositeSign': [
            ['Wjets_emu_BaryonMisMatched_OS', proxies.GenericChannel(region="BaryonMisMatched", name="BaryonMisMatchedOS")],
            ['SPG_BaryonMisMatched', proxies.GenericChannel(name="BaryonMisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS"])],
            # ['SPG_BaryonMisMatched', proxies.GenericChannel(name="BaryonMisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS", "-inclusive_Dplus_SS"])],  # noqa: E501
        ],
        'BaryonMisMatched_SameSign': [
            ['Wjets_emu_BaryonMisMatched_SS', proxies.GenericChannel(region="BaryonMisMatched", name="BaryonMisMatchedSS", os_minus_ss_fit_configuration=True, ss_only=True)],
            ['SPG_BaryonMisMatched', proxies.GenericChannel(name="BaryonMisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_SS"])],
            # ['SPG_BaryonMisMatched', proxies.GenericChannel(name="BaryonMisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS", "-inclusive_Dplus_SS"])],  # noqa: E501
        ],
        'CharmMisMatched': [
            ['Wjets_emu_CharmMisMatched_OS-SS', proxies.MatchedCharm(os_minus_ss_fit_configuration=True, name="CharmMisMatched_OS-SS")],
            # ['Wjets_emu_CharmMisMatched_OS', proxies.MatchedCharm(name="CharmMisMatched_OS")],
            # ['Wjets_emu_CharmMisMatched_SS', proxies.MatchedCharm(os_minus_ss_fit_configuration=True, ss_only=True, name="CharmMisMatched_SS")],
            ['SPG_CharmMisMatched', proxies.GenericChannel(name="CharmMisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS", "-inclusive_Dplus_SS"])],
        ],
        'CharmMisMatched_OppositeSign': [
            ['Wjets_emu_CharmMisMatched_OS', proxies.MatchedCharm(name="CharmMisMatched_OS")],
            ['SPG_CharmMisMatched', proxies.GenericChannel(name="CharmMisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_OS"])],
        ],
        'CharmMisMatched_SameSign': [
            ['Wjets_emu_CharmMisMatched_SS', proxies.MatchedCharm(os_minus_ss_fit_configuration=True, ss_only=True, name="CharmMisMatched_SS")],
            ['SPG_CharmMisMatched', proxies.GenericChannel(name="CharmMisMatchedInclusiveSPG", regions_override=["inclusive_Dplus_SS"])],
        ],
        # 'NoCharmBkg': [
        #     ['Wjets_Rest_OS-SS', proxies.NoMatchBackground(os_minus_ss_fit_configuration=True, add_no_truth_match=True)],
        #     ['Wjets_Rest_Loose_OS-SS', proxies.NoMatchBackground(os_minus_ss_fit_configuration=True, add_no_truth_match=True, loose_sr=True, regions_override=[
        #                                                              "el_minus_SR_0tag_Dplus_OS", "el_plus_SR_0tag_Dplus_OS", "mu_minus_SR_0tag_Dplus_OS", "mu_plus_SR_0tag_Dplus_OS"])],  # noqa: E501
        # ],
        # 'MisMatchBkg': [
        #     ['Wjets_emu_MisMatched_OS-SS', proxies.GenericChannel(name="MisMatched", os_minus_ss_fit_configuration=True, region=["MisMatched", "MatchedNoFid"])],  # noqa: E501
        #     ['Wjets_emu_MisMatched_Loose_OS-SS', proxies.GenericChannel(name="MisMatchedLoose", os_minus_ss_fit_configuration=True, region=["MisMatched", "MatchedNoFid"], regions_override=[  # noqa: E501
        #                                                                 "el_minus_SR_0tag_Dplus_OS", "el_plus_SR_0tag_Dplus_OS", "mu_minus_SR_0tag_Dplus_OS", "mu_plus_SR_0tag_Dplus_OS"])],  # noqa: E501
        # ],
    }

    def get(self):
        return self.samples


class SPGSamples:

    samples = {
        'Matched': [
            ['SPG_Matched', proxies.GenericChannel(name="MatchedSPG", regions_override=["inclusive_Dplus_OS_Matched", "-inclusive_Dplus_SS_Matched"])],
        ],
        '411MisMatched': [
            ['SPG_411MisMatched', proxies.GenericChannel(name="411MisMatchedSPG", regions_override=[
                                                         "inclusive_Dplus_OS_411MisMatched", "-inclusive_Dplus_SS_411MisMatched"])],
        ],
        'CharmMisMatched': [
            ['SPG_CharmMisMatched', proxies.GenericChannel(name="CharmMisMatchedSPG", regions_override=["inclusive_Dplus_OS", "-inclusive_Dplus_SS"])],
        ],
        'NoCharmBkg': [
            ['Wjets_Rest_Loose_OS-SS', proxies.NoMatchBackground(loose_sr=True, regions_override=[
                "el_minus_SR_0tag_Dplus_OS", "el_plus_SR_0tag_Dplus_OS", "mu_minus_SR_0tag_Dplus_OS", "mu_plus_SR_0tag_Dplus_OS",
                "-el_minus_SR_0tag_Dplus_SS", "-el_plus_SR_0tag_Dplus_SS", "-mu_minus_SR_0tag_Dplus_SS", "-mu_plus_SR_0tag_Dplus_SS"])],
        ],
        'MisMatchBkg': [
            ['Wjets_emu_MisMatched_Loose_OS-SS', proxies.GenericChannel(name="MisMatchedLoose", region=["MisMatched", "MatchedNoFid"], regions_override=[
                "el_minus_SR_0tag_Dplus_OS", "el_plus_SR_0tag_Dplus_OS", "mu_minus_SR_0tag_Dplus_OS", "mu_plus_SR_0tag_Dplus_OS",
                "-el_minus_SR_0tag_Dplus_SS", "-el_plus_SR_0tag_Dplus_SS", "-mu_minus_SR_0tag_Dplus_SS", "-mu_plus_SR_0tag_Dplus_SS"])],
        ],
    }

    def get(self):
        return self.samples


class WDComparisonSamples:

    samples = [
        ['MG_Wjets_emu_Matched', proxies.Matched()],
        ['Powheg_Wjets_emu_Matched', proxies.Matched()],
        ['Sherpa_Wjets_emu_Matched', proxies.Matched()],
        # ['STDM13_Wjets_emu_Matched', proxies.Matched()],
        # ['STDM13_Wjets_emu_411MisMatched', proxies.MisMatched(pdgId="411")],
        # ['FTAG4_Wjets_emu_Matched', proxies.Matched()],
        # ['FTAG4_Wjets_emu_411MisMatched', proxies.MisMatched(pdgId="411")],
    ]

    def get(self):
        return self.samples


class WDFlavourComparison:

    samples = [
        ['MG_Wjets_light_Rest', proxies.Rest()],
        ['MG_Wjets_cjets_Rest', proxies.Rest()],
        ['MG_Wjets_bjets_Rest', proxies.Rest()],
    ]

    def get(self):
        return self.samples


class WDBackgroundComparison:

    samples = [
        ['Wjets_bujets_emu_Rest_SR', proxies.Rest()],
        ['Wjets_bujets_emu_Rest_LooseSR', proxies.Rest(loose_sr=True, name="Loose")],
        ['Wjets_cjets_emu_Rest_SR', proxies.Rest()],
        ['Wjets_cjets_emu_Rest_LooseSR', proxies.Rest(loose_sr=True, name="Loose")],
    ]

    def get(self):
        return self.samples


class ChannelGenerator:

    def __init__(self, config, samples, signs, years, leptons, charges,
                 btags, ptbins=[""], sample_config="", decay_mode="", process_string="",
                 force_positive=False, replacement_samples={}, os_minus_ss_fit_configuration=False, make_plots=True, save_to_file=True):
        self.config = config
        self.decay_mode = decay_mode
        self.signs = signs
        self.years = years
        self.leptons = leptons
        self.charges = charges
        self.btags = btags
        self.ptbins = ptbins
        self.process_string = process_string
        self.sample_config = sample_config
        self.force_positive = force_positive
        self.replacement_samples = replacement_samples
        self.template = samples
        self.make_plots = make_plots
        self.save_to_file = save_to_file
        self.os_minus_ss_fit_configuration = os_minus_ss_fit_configuration
        if type(samples.get()) == dict:
            self.samples = samples.get()
        else:
            self.samples = {'': samples.get()}
        if type(self.years) == list and len(self.years) == 0:
            self.years = ['']

    def get_config(self):
        return self.config

    def make_channel(self, lumi, sign='', year='', lepton='', charge='', btag='', extra_rebin=1, os_only=False):
        for ptbin in self.ptbins:
            for suffix, samples in self.samples.items():
                channel_name = self.generate_channel_name(sign=sign, year=year, lepton=lepton, charge=charge, btag=btag, suffix=suffix)
                if ptbin:
                    channel_name += f"_{ptbin}"
                regions = self.generate_channel_regions(sign=sign, year=year, lepton=lepton, charge=charge, btag=btag)
                if ptbin and ptbin != "inc":
                    regions_modified = [f"{x}_{ptbin}" for x in regions]
                elif ptbin == "inc":
                    regions_modified = [f"{x}_{pt}" for x in regions for pt in self.ptbins]
                channel = {
                    'regions': regions_modified if ptbin else regions,
                    'label': self.generate_channel_labels(sign=sign, lepton=lepton, charge=charge, btag=btag, extra=ptbin),
                    'lumi': '+'.join(lumi),
                    'extra_rebin': extra_rebin,
                    'force_positive': self.force_positive,
                    'save_to_file': self.save_to_file,
                    'samples': [],
                }

                # SPG samples
                if self.os_minus_ss_fit_configuration:
                    if sign == 'OS' and btag != '1tag':
                        channel['replacement_samples'] = {k: v for k, v in self.replacement_samples.items()}
                else:
                    if sign == '' and btag != '1tag':
                        channel['replacement_samples'] = {k: v for k, v in self.replacement_samples.items()}

                # pt bins for SPG samples
                if 'replacement_samples' in channel:
                    for k, v in channel['replacement_samples'].items():
                        if ptbin:
                            channel['replacement_samples'][k] = f"{v}_{ptbin}"
                        else:
                            channel['replacement_samples'][k] = f"{v}"

                # add to config
                self.config['channels'][channel_name] = channel

                # make plots
                if not self.make_plots:
                    channel['make_plots'] = False

                # scale factors
                if hasattr(self.template, "scale_factors"):
                    channel["scale_factors"] = self.template.scale_factors

                # print out
                print(f'Generated channel {channel_name} with {len(regions)} regions')
                for reg in regions:
                    print(f'  {reg}')

                # add samples
                for sample in samples:
                    if len(sample) > 1 and sample[1].os_minus_ss_fit_configuration and 'SS' in channel_name:
                        continue
                    if len(sample) == 1:
                        channel['samples'] += [sample[0]]
                    else:
                        # generate proxy channel
                        proxy = sample[1]
                        proxy_channel_name = channel_name + "_" + proxy.get_name()
                        if ptbin and ptbin != "inc":
                            proxy_regions_modified = [f"{x}_{ptbin}" for x in proxy.get_regions(regions)]
                        elif ptbin == "inc":
                            proxy_regions_modified = [f"{x}_{pt}" for x in proxy.get_regions(regions) for pt in self.ptbins]
                        proxy_channel = {
                            'make_plots': False,
                            'save_to_file': False,
                            'regions': proxy_regions_modified if ptbin else proxy.get_regions(regions),
                        }
                        self.config['channels'][proxy_channel_name] = proxy_channel
                        channel['samples'] += [f'{sample[0]} | {proxy_channel_name}']

    def generate_channel_name(self, sign='', year='', lepton='', charge='', btag='', suffix=''):
        if sign:
            sign = f'{sign}_'
        else:
            sign = 'OS-SS_'
        if suffix:
            suffix = f'_{suffix}'
        btag_massaged = [btag] if type(btag) != list else btag
        name = f'{sign}{format(year)}{format(lepton)}{format(charge)}{format(btag_massaged[0])}{self.decay_mode}{suffix}'
        return name

    def generate_channel_regions(self, sign='', year='', lepton='', charge='', btag=''):
        regions = []
        if not btag:
            btag = flatten(self.btags)
        btag_massaged = [btag] if type(btag) != list else btag
        if sign:
            regions = [f'{format(y)}{format(lep)}{format(c)}SR_{format(b)}{self.decay_mode}_{sign}'
                       for y in (self.years if not year else [year])
                       for lep in (self.leptons if not lepton else [lepton])
                       for c in (self.charges if not charge else [charge])
                       for b in (self.btags if not btag else btag_massaged)]
        else:
            regions = [f'{format(y)}{format(lep)}{format(c)}SR_{format(b)}{self.decay_mode}_OS'
                       for y in (self.years if not year else [year])
                       for lep in (self.leptons if not lepton else [lepton])
                       for c in (self.charges if not charge else [charge])
                       for b in (self.btags if not btag else btag_massaged)] + \
                [f'-{format(y)}{format(lep)}{format(c)}SR_{format(b)}{self.decay_mode}_SS'
                 for y in (self.years if not year else [year])
                 for lep in (self.leptons if not lepton else [lepton])
                 for c in (self.charges if not charge else [charge])
                 for b in (self.btags if not btag else btag_massaged)]
        return regions

    def generate_channel_labels(self, sign='', lepton='', charge='', btag='', extra=''):
        labels = [f'{self.process_string}, {sign if sign else "OS-SS"}']
        row2 = ''
        if charge:
            charge = f' {charge}'
        if lepton:
            row2 = f'{lepton}{charge} channel'
        else:
            row2 = 'inclusive channel'
        if btag:
            btag_massaged = [btag] if type(btag) != list else btag
            row2 += f', {"+".join(btag_massaged)}'
        labels += [row2]
        if self.sample_config:
            labels += [label_dict[self.sample_config]]
        if extra:
            labels += [extra]
        return labels
