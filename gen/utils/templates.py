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
        if self.sample_config not in ['truth_comparison', 'wplusd_comparison', 'spg_comparison']:
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

    def __init__(self, os_minus_ss_fit_configuration=False, loose_sr=False):
        self.os_minus_ss_fit_configuration = os_minus_ss_fit_configuration
        self.loose_sr = loose_sr
        self.samples = []
        if self.os_minus_ss_fit_configuration:
            self.samples += [['MockMC', proxies.MockMC()]]
        self.samples += [
            ['Wjets_emu_Matched', proxies.Matched(os_minus_ss_fit_configuration=self.os_minus_ss_fit_configuration)],
            ['Wjets_emu_411MisMatched', proxies.MisMatched(
                os_minus_ss_fit_configuration=self.os_minus_ss_fit_configuration, loose_sr=self.loose_sr, pdgId="411")],
            ['Wjets_emu_Charm', proxies.MatchedCharm(os_minus_ss_fit_configuration=self.os_minus_ss_fit_configuration, loose_sr=self.loose_sr)],
            ['Wjets_emu_MisMatched', proxies.GenericChannel(
                name="MisMatched", os_minus_ss_fit_configuration=self.os_minus_ss_fit_configuration, region=["MisMatched", "MatchedNoFid"])],
            ['Wjets_emu_Rest', proxies.NoMatchBackground(os_minus_ss_fit_configuration=self.os_minus_ss_fit_configuration, loose_sr=self.loose_sr)],
            ['Top_Matched', proxies.Matched(os_minus_ss_fit_configuration=self.os_minus_ss_fit_configuration)],
            ['Top_Rest', proxies.Rest(os_minus_ss_fit_configuration=self.os_minus_ss_fit_configuration, include_other=True, name="WithOther")],
            ['DibosonZjetsTau', proxies.PlainChannel(os_minus_ss_fit_configuration=self.os_minus_ss_fit_configuration)],
            ['Multijet_MatrixMethod', proxies.MatrixMethod(os_minus_ss_fit_configuration=self.os_minus_ss_fit_configuration)],
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
            ['SPG_DPlus_Matched', proxies.GenericChannel(name="MatchedSPG", regions_override=["inclusive_Dplus_per_D_Matched", "-inclusive_Dplus_per_antiD_Matched"])],
        ],
        '411MisMatched': [
            ['Wjets_emu_411MisMatched_OS-SS', proxies.GenericChannel(region="411MisMatched", name="411MisMatchedOS-SS", os_minus_ss_fit_configuration=True)],
            ['SPG_DPlus_411MisMatched', proxies.GenericChannel(name="411MisMatchedSPG", regions_override=["inclusive_Dplus_per_D_411MisMatched", "-inclusive_Dplus_per_antiD_411MisMatched"])],
        ],
        '431MisMatched': [
            ['Wjets_emu_431MisMatched_OS-SS', proxies.GenericChannel(region="431MisMatched", name="431MisMatchedOS-SS", os_minus_ss_fit_configuration=True)],
            ['SPG_DPlus_431MisMatched', proxies.GenericChannel(name="431MisMatchedSPG", regions_override=["inclusive_Dplus_per_D", "-inclusive_Dplus_per_antiD"])],
        ],
        '421MisMatched': [
            ['Wjets_emu_421MisMatched_OS-SS', proxies.GenericChannel(region=["421MisMatched", "413MisMatched"], name="421MisMatchedOS-SS", os_minus_ss_fit_configuration=True)],
            ['SPG_DPlus_421MisMatched', proxies.GenericChannel(name="421MisMatchedSPG", regions_override=["inclusive_Dplus_per_D", "-inclusive_Dplus_per_antiD"])],
        ],
        'BaryonMisMatched': [
            ['Wjets_emu_BaryonMisMatched_OS-SS', proxies.GenericChannel(region="BaryonMisMatched", name="BaryonMisMatchedOS-SS", os_minus_ss_fit_configuration=True)],
            ['SPG_DPlus_BaryonMisMatched', proxies.GenericChannel(name="BaryonMisMatchedSPG", regions_override=["inclusive_Dplus_per_D", "-inclusive_Dplus_per_antiD"])],
        ],
        'CharmMisMatched': [
            ['Wjets_emu_CharmMisMatched_OS-SS', proxies.MatchedCharm(os_minus_ss_fit_configuration=True)],
            ['SPG_DPlus_CharmMisMatched', proxies.GenericChannel(name="CharmMisMatchedSPG", regions_override=["inclusive_Dplus_per_D", "-inclusive_Dplus_per_antiD"])],
        ],
    }

    scale_factors = {
        'SPG_Ds': 0.2991,
        'SPG_D0': 0.1500,
        'SPG_Baryon': 0.2176,
    }

    def get(self):
        return self.samples


class SPGSamples:

    samples = {
        'Matched': [
            ['SPG_DPlus_Matched', proxies.GenericChannel(name="MatchedSPG", regions_override=["inclusive_Dplus_per_D_Matched", "-inclusive_Dplus_per_antiD_Matched"])],
        ],
        '411MisMatched': [
            ['SPG_DPlus_411MisMatched', proxies.GenericChannel(name="411MisMatchedSPG", regions_override=["inclusive_Dplus_per_D_411MisMatched", "-inclusive_Dplus_per_antiD_411MisMatched"])],
        ],
        'CharmMisMatched': [
            ['SPG_DPlus_CharmMisMatched', proxies.GenericChannel(name="CharmMisMatchedSPG", regions_override=["inclusive_Dplus_per_D", "-inclusive_Dplus_per_antiD"])],
        ],
    }

    scale_factors = {
        'SPG_Ds': 0.2991,
        'SPG_D0': 0.1500,
        'SPG_Baryon': 0.2176,
    }

    def get(self):
        return self.samples

class WDComparisonSamples:

    samples = [
        ['MG_Wjets_emu_Rest_el_minus', proxies.Rest(allowed_regions=["el_minus"], name="el_minus")],
        ['MG_Wjets_emu_Rest_el_plus', proxies.Rest(allowed_regions=["el_plus"], name="el_plus")],
        ['MG_Wjets_emu_Rest_mu_minus', proxies.Rest(allowed_regions=["mu_minus"], name="mu_minus")],
        ['MG_Wjets_emu_Rest_mu_plus', proxies.Rest(allowed_regions=["mu_plus"], name="mu_plus")],
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
                 btags, sample_config="", decay_mode="", process_string="",
                 force_positive = False, replacement_samples = {}, make_plots=True):
        self.config=config
        self.decay_mode=decay_mode
        self.signs=signs
        self.years=years
        self.leptons=leptons
        self.charges=charges
        self.btags=btags
        self.process_string=process_string
        self.sample_config=sample_config
        self.force_positive=force_positive
        self.replacement_samples=replacement_samples
        self.template=samples
        self.make_plots=make_plots
        if type(samples.get()) == dict:
            self.samples=samples.get()
        else:
            self.samples={'': samples.get()}
        if type(self.years) == list and len(self.years) == 0:
            self.years=['']

    def get_config(self):
        return self.config

    def make_channel(self, lumi, sign = '', year = '', lepton = '', charge = '', btag = '', extra_rebin = 1, os_only = False):
        for suffix, samples in self.samples.items():
            channel_name=self.generate_channel_name(sign = sign, year = year, lepton = lepton, charge = charge, btag = btag, suffix = suffix)
            regions=self.generate_channel_regions(sign = sign, year = year, lepton = lepton, charge = charge, btag = btag)
            channel={
                'regions': regions,
                'label': self.generate_channel_labels(sign=sign, lepton=lepton, charge=charge, btag=btag),
                'lumi': '+'.join(lumi),
                'extra_rebin': extra_rebin,
                'force_positive': self.force_positive,
                'save_to_file': True,
                'samples': [],
            }
            if sign == 'OS':
                channel['replacement_samples'] = self.replacement_samples
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
                if 'MockMC' not in sample[0] and os_only and 'SS' in channel_name:
                    continue
                if len(sample) == 1:
                    channel['samples'] += [sample[0]]
                else:
                    # generate proxy channel
                    proxy = sample[1]
                    proxy_channel_name = channel_name + "_" + proxy.get_name()
                    proxy_channel = {
                        'make_plots': False,
                        'regions': proxy.get_regions(regions),
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

    def generate_channel_labels(self, sign='', lepton='', charge='', btag=''):
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
        return labels
