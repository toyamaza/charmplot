#!/usr/bin/env python
import abc
import gen.utils.proxies as proxies


def flatten(xs):
    result = []
    if isinstance(xs, (list, tuple)):
        for x in xs:
            result.extend(flatten(x))
    else:
        result.append(xs)
    return result


def format(string):
    if string:
        return f'{string}_'
    else:
        return string


class DataMCConfig:

    def __init__(self, variables, systematics=[]):
        self.variables = variables
        self.systematics = systematics

    def to_dict(self):
        out = {
            'variablesConf': self.variables,
            'channels': {},
        }
        if self.systematics:
            out['systematics'] = self.systematics
        return out


class ChannelTemplate:
    labels = []
    samplesConf = None
    data = None

    @property
    @abc.abstractmethod
    def get(self):
        pass


class WDTruthSamples(ChannelTemplate):

    # base class
    samplesConf = "madgraph_truth"
    data = "Data"

    # new objects
    truthSlices = [
        "Matched_truth_pt_bin1",
        "Matched_truth_pt_bin2",
        "Matched_truth_pt_bin3",
        "Matched_truth_pt_bin4",
        "Matched_truth_pt_bin5",
    ]

    def __init__(self, fitType="", MockMC=True, decayMode="Dplus", truthDiffBins=False, splitSignalSamples=False):
        self.os_ss_sub = fitType == "OS-SS"
        self.MockMC = MockMC
        self.decayMode = decayMode
        self.truthDiffBins = truthDiffBins
        self.splitSignalSamples = splitSignalSamples
        self.samples = []

        if self.decayMode == "Dstar":
            self.samplesConf = "madgraph_truth_dstar"

        # MockMC at the top in case of OS-SS plots
        if self.os_ss_sub and self.MockMC:
            self.samples += [['MockMC', proxies.MockMC(subtract_mj=False)]]

        # signal sample
        if self.truthDiffBins and self.splitSignalSamples:
            self.samples += [[f'Wjets_emu_{slice}', proxies.Matched(os_ss_sub=self.os_ss_sub, ptbin=i + 1)] for i, slice in enumerate(self.truthSlices)]
        else:
            if self.splitSignalSamples:
                self.samples += [['Wjets_emu_Matched', proxies.GenericChannel(region=self.truthSlices, name="Matched", os_ss_sub=self.os_ss_sub)]]
            else:
                self.samples += [['Wjets_emu_Matched', proxies.Matched(os_ss_sub=self.os_ss_sub)]]

        # backgrdound from other than the signal decay modes
        self.samples += [
            ['Wjets_emu_Charm', proxies.GenericChannel(region=["411MisMatched", "413MisMatched", "421MisMatched", "431MisMatched", "BaryonMisMatched", "Wjets_emu_Charm"], name="MatchedCharm")],
            ['Wjets_emu_MisMatched', proxies.GenericChannel(name="MisMatched", os_ss_sub=self.os_ss_sub, region=["MisMatched", "MatchedNoFid"])],
            ['Wjets_emu_Rest', proxies.GenericChannel(name="Rest", os_ss_sub=self.os_ss_sub, region=["Other", "HardMisMatched"])],
            ['Top', proxies.PlainChannel(os_ss_sub=self.os_ss_sub)],
            ['DibosonVjetsTau', proxies.PlainChannel(os_ss_sub=self.os_ss_sub)],
        ]

        # MultiJet
        self.samples += [['Multijet_MatrixMethod', proxies.MatrixMethod(os_ss_sub=self.os_ss_sub, fake_factor=False)]]

        # MockMC at the bottom
        if not self.os_ss_sub and self.MockMC:
            self.samples += [['MockMC_minus_MC', proxies.MockMC(subtract_mj=True)]]
            self.samples += [['Offset', proxies.PlainChannel()]]

    def get(self):
        return self.samples


class WDFlavourSamples(ChannelTemplate):

    samplesConf = "madgraph"
    data = "Data"

    samples = [
        ['Wjets_cjets_emu'],
        ['Wjets_bjets_emu'],
        ['Wjets_light_emu'],
        ['Top'],
        ['DibosonVjetsTau'],
        ['Multijet_MatrixMethod', proxies.MatrixMethod()]
    ]

    def get(self):
        return self.samples


class WDFakeTrackSamples(ChannelTemplate):

    samplesConf = "madgraph"

    samples = [
        ['Zero_Fake', proxies.FakeTrack0()],
        ['One_Fake', proxies.FakeTrack1()],
        ['Two_Fake', proxies.FakeTrack2()],
        ['Three_Fake', proxies.FakeTrack3()]
    ]

    def get(self):
        return self.samples


class SPGComparison(ChannelTemplate):

    # base class
    samplesConf = "spg_comparison"

    # new objects
    truthSlices = [
        "Matched_truth_pt_bin1",
        "Matched_truth_pt_bin2",
        "Matched_truth_pt_bin3",
        "Matched_truth_pt_bin4",
        "Matched_truth_pt_bin5",
    ]

    def __init__(self, truthDiffBins=False, splitSignalSamples=False, decay_mode="Dplus", signal_only=False, background_only=False):
        self.truthDiffBins = truthDiffBins
        self.splitSignalSamples = splitSignalSamples
        self.decay_mode = decay_mode
        self.signal_only = signal_only
        self.background_only = background_only

        if self.decay_mode == "Dstar":
            self.samplesConf = "spg_comparison_dstar"

        self.samples = {}

        # signal samples
        if not self.background_only:
            if self.splitSignalSamples:
                self.samples.update(
                    {
                        'Matched': [
                            ['Wjets_emu_Matched', proxies.GenericChannel(region=self.truthSlices, name="Matched")],
                            ['SPG_Matched', proxies.SPGChannel(name="SPG_Matched",
                                                               regions_OS=["inclusive_" + self.decay_mode + f"_OS_{slice}" for slice in self.truthSlices],
                                                               regions_SS=["inclusive_" + self.decay_mode + f"_SS_{slice}" for slice in self.truthSlices],
                                                               always_OS=True)]]
                    }
                )
            else:
                self.samples.update(
                    {
                        'Matched': [
                            ['Wjets_emu_Matched', proxies.Matched()],
                            ['SPG_Matched', proxies.SPGChannel(name="SPG_Matched",
                                                               regions_OS=["inclusive_" + self.decay_mode + "_OS_Matched"],
                                                               regions_SS=["inclusive_" + self.decay_mode + "_SS_Matched"],
                                                               always_OS=True)]]
                    }
                )

            # signal samples in truth differential bins
            if self.truthDiffBins and self.splitSignalSamples:
                self.samples.update(
                    {
                        slice: [
                            [f'Wjets_emu_{slice}', proxies.GenericChannel(region=slice, name=slice)],
                            [f'SPG_{slice}', proxies.SPGChannel(name=f"SPG_{slice}",
                                                                regions_OS=["inclusive_" + self.decay_mode + f"_OS_{slice}"],
                                                                regions_SS=["inclusive_" + self.decay_mode + f"_SS_{slice}"],
                                                                always_OS=True)],
                        ] for slice in self.truthSlices
                    })

        if self.decay_mode == "Dplus":
            if not self.signal_only:
                self.samples.update(
                    {
                        '411MisMatched': [
                            ['Wjets_emu_411MisMatched', proxies.GenericChannel(region="411MisMatched", name="411MisMatched")],
                            ['SPG_411MisMatched', proxies.SPGChannel(name="411MisMatchedInclusiveSPG", regions_OS=["inclusive_Dplus_OS_411MisMatched"], regions_SS=[
                                "inclusive_Dplus_SS_411MisMatched"], always_OS=True)],
                        ],
                        '421MisMatched': [
                            ['Wjets_emu_421MisMatched', proxies.GenericChannel(region=["421MisMatched", "413MisMatched"], name="421MisMatched")],
                            ['SPG_421MisMatched', proxies.SPGChannel(name="421MisMatchedInclusiveSPG", regions_OS=[
                                "inclusive_Dplus_OS"], regions_SS=["inclusive_Dplus_SS"])],
                        ],
                    }
                )
        elif self.decay_mode == "Dstar":
            if not self.signal_only:
                self.samples.update(
                    {
                        '413MisMatched': [
                            ['Wjets_emu_413MisMatched', proxies.GenericChannel(region="413MisMatched", name="413MisMatched")],
                            ['SPG_413MisMatched', proxies.SPGChannel(name="413MisMatchedInclusiveSPG", regions_OS=[
                                                                     "inclusive_Dstar_OS_413MisMatched"], regions_SS=["inclusive_Dstar_SS_413MisMatched"])],
                        ],
                        '411MisMatched': [
                            ['Wjets_emu_411MisMatched', proxies.GenericChannel(region="411MisMatched", name="411MisMatched")],
                            ['SPG_411MisMatched', proxies.SPGChannel(name="411MisMatchedInclusiveSPG", regions_OS=[
                                                                     "inclusive_Dstar_OS"], regions_SS=["inclusive_Dstar_SS"], always_OS=True)],
                        ],
                        '421MisMatched': [
                            ['Wjets_emu_421MisMatched', proxies.GenericChannel(region="421MisMatched", name="421MisMatched")],
                            ['SPG_421MisMatched', proxies.SPGChannel(name="421MisMatchedInclusiveSPG", regions_OS=[
                                                                     "inclusive_Dstar_OS"], regions_SS=["inclusive_Dstar_SS"])],
                        ],
                    }
                )

        # backgrounds
        if not self.signal_only:
            self.samples.update(
                {
                    '431MisMatched': [
                        ['Wjets_emu_431MisMatched', proxies.GenericChannel(region="431MisMatched", name="431MisMatched")],
                        ['SPG_431MisMatched', proxies.SPGChannel(name="431MisMatchedInclusiveSPG", regions_OS=[
                            "inclusive_" + self.decay_mode + "_OS"], regions_SS=["inclusive_" + self.decay_mode + "_SS"])],
                    ],
                    'BaryonMisMatched': [
                        ['Wjets_emu_BaryonMisMatched', proxies.GenericChannel(region="BaryonMisMatched", name="BaryonMisMatched")],
                        ['SPG_BaryonMisMatched', proxies.SPGChannel(name="BaryonMisMatchedInclusiveSPG", regions_OS=[
                            "inclusive_" + self.decay_mode + "_OS"], regions_SS=["inclusive_" + self.decay_mode + "_SS"])],
                    ],
                    'CharmMisMatched': [
                        ['Wjets_emu_CharmMisMatched', proxies.MatchedCharm(name="CharmMisMatched")],
                        ['SPG_CharmMisMatched', proxies.SPGChannel(name="CharmMisMatchedInclusiveSPG", regions_OS=[
                            "inclusive_" + self.decay_mode + "_OS"], regions_SS=["inclusive_" + self.decay_mode + "_SS"])],
                    ],
                })

    def get(self):
        return self.samples


class BKGComparison(ChannelTemplate):

    samplesConf = "bkg_comparison"

    samples = {
        'Wjets_emu_Rest': [
            ['Wjets_emu_Rest', proxies.NoMatchBackground()],
            ['Wjets_emu_Rest_PostProc', proxies.GenericChannel(name="Wjets_emu_Rest", regions_OS=["Wjets_emu_Rest_OS"], regions_SS=["Wjets_emu_Rest_SS"])],
        ],
        'Wjets_emu_MisMatched': [
            ['Wjets_emu_MisMatched', proxies.GenericChannel(name="MisMatched", region=["MisMatched", "MatchedNoFid"])],
            ['Wjets_emu_MisMatched_PostProc', proxies.GenericChannel(name="Wjets_emu_MisMatched", regions_OS=[
                                                                     "Wjets_emu_MisMatched_OS"], regions_SS=["Wjets_emu_MisMatched_SS"])],
        ],
        # 'Wjets_emu_Rest': [
        #     ['Wjets_emu_Rest', proxies.NoMatchBackground()],
        #     ['Sherpa_Wjets_emu_Rest', proxies.NoMatchBackground()],
        #     ['Wjets_emu_Rest_PostProc', proxies.GenericChannel(name="Wjets_emu_Rest", regions_OS=["Wjets_emu_Rest_OS"], regions_SS=["Wjets_emu_Rest_SS"])],
        #     ['Sherpa_Wjets_emu_Rest_PostProc', proxies.GenericChannel(name="Wjets_emu_Rest", regions_OS=["Wjets_emu_Rest_OS"], regions_SS=["Wjets_emu_Rest_SS"])],
        # ],
        # 'DibosonVjetsTau': [
        #     ['DibosonVjetsTau', proxies.PlainChannel()],
        #     ['DibosonVjetsTau_PostProc', proxies.GenericChannel(name="DibosonVjetsTau", regions_OS=["Other_OS"], regions_SS=["Other_SS"])],
        # ],
        # 'MultiJet': [
        #     ['Multijet_MatrixMethod', proxies.MatrixMethod()],
        #     ['MultiJet_PostProc', proxies.GenericChannel(name="MultiJet", regions_OS=["MultiJet_OS"], regions_SS=["MultiJet_SS"])],
        # ],
    }

    def get(self):
        return self.samples


class SPGSysComparison(ChannelTemplate):

    samplesConf = "madgraph_truth"

    samples = {
        'MatchedCharm': [
            ['Wjets_emu_Charm', proxies.GenericChannel(region=["411MisMatched", "413MisMatched", "421MisMatched", "431MisMatched", "BaryonMisMatched"], name="MatchedCharm")],
            # ['SPG_CharmMisMatched_sys', proxies.GenericChannel(region=["Wjets_emu_Charm"], name="SPGMatchedCharm")],
            ['SPG_CharmMisMatched_sys_Dplus', proxies.GenericChannel(region=["Wjets_emu_Charm"], name="SPGMatchedCharmDplus")],
            ['SPG_CharmMisMatched_sys_Dzero', proxies.GenericChannel(region=["Wjets_emu_Charm"], name="SPGMatchedCharmDzero")],
            ['SPG_CharmMisMatched_sys_Dsubs', proxies.GenericChannel(region=["Wjets_emu_Charm"], name="SPGMatchedCharmDsubs")],
        ],
    }

    def get(self):
        return self.samples


class ReplacementSamples(ChannelTemplate):

    truthSlices = [
        "Matched_truth_pt_bin1",
        "Matched_truth_pt_bin2",
        "Matched_truth_pt_bin3",
        "Matched_truth_pt_bin4",
        "Matched_truth_pt_bin5",
    ]

    def __init__(self, truthDiffBins=False, splitSignalSamples=False, decay_mode="Dplus"):
        self.truthDiffBins = truthDiffBins
        self.splitSignalSamples = splitSignalSamples
        self.decay_mode = decay_mode

        # signal samples
        if self.splitSignalSamples:
            self.samples = {
                'Matched': [['SPG_Matched', proxies.SPGChannel(name="SPG_Matched",
                                                               regions_OS=[f"inclusive_{self.decay_mode}_OS_{slice}" for slice in self.truthSlices],
                                                               regions_SS=[f"inclusive_{self.decay_mode}_SS_{slice}" for slice in self.truthSlices],
                                                               always_OS=True)]]
            }
        else:
            self.samples = {
                'Matched': [['SPG_Matched', proxies.SPGChannel(name="SPG_Matched",
                                                               regions_OS=[f"inclusive_{self.decay_mode}_OS_Matched"],
                                                               regions_SS=[f"inclusive_{self.decay_mode}_SS_Matched"],
                                                               always_OS=True)]]
            }

        # signal samples in truth differential bins
        if self.truthDiffBins and self.splitSignalSamples:
            self.samples.update(
                {
                    slice: [[f'SPG_{slice}', proxies.SPGChannel(name=f"SPG_{slice}",
                                                                regions_OS=[f"inclusive_{self.decay_mode}_OS_{slice}"],
                                                                regions_SS=[f"inclusive_{self.decay_mode}_SS_{slice}"],
                                                                always_OS=True)],
                            ] for slice in self.truthSlices
                })

        self.samples.update({
            'CharmMisMatched': [
                ['SPG_CharmMisMatched', proxies.SPGChannel(name="CharmMisMatchedInclusiveSPG", regions_OS=[
                    f"inclusive_{self.decay_mode}_OS"], regions_SS=[f"inclusive_{self.decay_mode}_SS"])],
            ],
            'Wjets_emu_MisMatched': [
                ['Wjets_emu_MisMatched_PostProc', proxies.GenericChannel(name="Wjets_emu_MisMatched", regions_OS=[
                    "Wjets_emu_MisMatched_OS"], regions_SS=["Wjets_emu_MisMatched_SS"])],
            ],
            'Wjets_emu_Rest': [
                ['Wjets_emu_Rest_PostProc', proxies.GenericChannel(name="Wjets_emu_Rest", regions_OS=[
                    "Wjets_emu_Rest_OS"], regions_SS=["Wjets_emu_Rest_SS"])],
                ['Sherpa_Wjets_emu_Rest_PostProc', proxies.GenericChannel(name="Wjets_emu_Rest", regions_OS=[
                    "Wjets_emu_Rest_OS"], regions_SS=["Wjets_emu_Rest_SS"])],
            ],
            # 'Wjets_emu_Bkg': [
            #     ['Wjets_emu_Bkg_PostProc', proxies.GenericChannel(name="Wjets_emu_Bkg", regions_OS=[
            #         "Wjets_emu_Bkg_OS"], regions_SS=["Wjets_emu_Bkg_SS"])],
            #     ['Sherpa_Wjets_emu_Bkg_PostProc', proxies.GenericChannel(name="Wjets_emu_Bkg", regions_OS=[
            #         "Wjets_emu_Bkg_OS"], regions_SS=["Wjets_emu_Bkg_SS"])],
            # ],
            # 'DibosonVjetsTau': [
            #     ['DibosonVjetsTau_PostProc', proxies.GenericChannel(name="DibosonVjetsTau", regions_OS=["Other_OS"], regions_SS=["Other_SS"])],
            # ],
        })

    def get(self):
        return self.samples


class WDComparisonSamples(ChannelTemplate):

    samplesConf = "wplusd_comparison"

    samples = [
        ['MG_Wjets_emu_Matched', proxies.Matched()],
        # ['Powheg_Wjets_emu_Matched', proxies.Matched()],
        ['Sherpa_Wjets_emu_Matched', proxies.Matched()],
        # ['MGFxFx_Wjets_emu_Matched', proxies.Matched()],
        # ['Sherpa2211_Wjets_emu_Matched', proxies.Matched()],
    ]

    def get(self):
        return self.samples


class WjetsSherpaSys(ChannelTemplate):

    samplesConf = "madgraph_truth"

    samples = [
        ['Wjets_emu_Rest', proxies.GenericChannel(name="Rest", region=["Other", "HardMisMatched"])],
        ['Sherpa_Wjets_emu_Rest', proxies.GenericChannel(name="Rest", region=["Other", "HardMisMatched"])],
        # ['Wjets_emu_Rest_Default', proxies.GenericChannel(name="Rest", region=["Other", "HardMisMatched"])],
        # ['Sherpa_Wjets_emu_Rest_Default', proxies.GenericChannel(name="Rest", region=["Other", "HardMisMatched"])],
    ]

    def get(self):
        return self.samples


class SignalComparison(ChannelTemplate):

    samplesConf = "madgraph_truth"

    samples = {}

    truthSlices = [
        "Matched_truth_pt_bin1",
        "Matched_truth_pt_bin2",
        "Matched_truth_pt_bin3",
        "Matched_truth_pt_bin4",
        "Matched_truth_pt_bin5",
    ]

    samples.update(
        {
            'Matched': [
                ['Wjets_emu_Matched_1tag', proxies.GenericChannel(region=truthSlices, name="Matched_1tag", force_1tag=True, loose_sr=True)],
                # ['Wjets_emu_Matched', proxies.GenericChannel(region=truthSlices, name="Matched")],
                ['Wjets_emu_Matched_Loose', proxies.GenericChannel(region=truthSlices, name="Matched_Loose", loose_sr=True)],
            ]
        }
    )

    def get(self):
        return self.samples


class ChannelGenerator:

    def __init__(self, config, samples, signs, years, leptons, charges,
                 btags, ptbins=[""], sample_config="", decay_mode="", process_string="",
                 force_positive=False, replacement_samples={}, os_ss_sub=False, make_plots=True, save_to_file=True):
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
        self.os_ss_sub = os_ss_sub
        if type(samples.get()) == dict:
            self.samples = samples.get()
        else:
            self.samples = {'': samples.get()}
        if type(self.years) == list and len(self.years) == 0:
            self.years = ['']

    def get_config(self):
        if self.template.samplesConf:
            self.config["samplesConf"] = self.template.samplesConf
        if self.template.data:
            self.config["data"] = self.template.data
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
                if btag == '0tag':
                    if channel_name.startswith("OS-SS"):
                        replacement_channel_sign = "OS-SS"
                    elif channel_name.startswith("OS"):
                        replacement_channel_sign = "OS"
                    elif channel_name.startswith("SS"):
                        replacement_channel_sign = "SS"
                    else:
                        raise Exception(f"Unrecognised charge in {channel_name}")
                    channel['replacement_samples'] = {k: v.replace("<charge>", replacement_channel_sign) for k, v in self.replacement_samples.items()}

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
                    if len(sample) > 1 and sample[1].os_ss_sub and 'SS' in channel_name:
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
        if lepton == '' and charge != '':
            lepton = 'lep'
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
        if lepton == "el":
            row2 = "e"
        elif lepton == "mu":
            row2 = "#mu"
        else:
            row2 = "W"
        if charge == "plus":
            row2 += "^{+}"
        elif charge == "minus":
            row2 += "^{-}"
        else:
            row2 += "^{#pm}"
        row2 += " channel"
        if btag:
            btag_massaged = [btag] if type(btag) != list else btag
            row2 += f', {"+".join(btag_massaged)}'
        labels += [row2]
        if self.template:
            labels += [x for x in self.template.labels]
        if extra:
            labels += [extra]
        return labels
