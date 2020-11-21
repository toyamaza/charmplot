#!/usr/bin/env python
import gen.utils.proxies as proxies

label_dict = {
    'madgraph_truth': 'MadGraph LO',
    'sherpa_truth': 'Sherpa 2.2.1',
    'sherpa2210_truth': 'Sherpa 2.2.10',
    'madgraph_fxfx_truth': 'MadGraph FxFx',
    'madgraph': 'MadGraph LO',
    'sherpa': 'Sherpa 2.2.1',
    'sherpa2210': 'Sherpa 2.2.10',
    'madgraph_fxfx': 'MadGraph FxFx',
}


def format(string):
    if string:
        return f'{string}_'
    else:
        return string


class DataMCConfig:

    def __init__(self, variables, samples):
        self.variables = variables
        self.samples = samples

    def to_dict(self):
        return {
            'variablesConf': self.variables,
            'samplesConf': self.samples,
            'data': 'Data',
            'channels': {},
        }


class WDTruthSamples:

    samples = [
        ['Wjets_emu_Matched', proxies.Matched()],
        ['Wjets_emu_Rest', proxies.Rest()],
        ['Top_Matched', proxies.Matched()],
        ['Top_Rest', proxies.Rest()],
        ['Wjets_emu_NoMatch', proxies.NoMatch()],
        ['Zjets_emu'],
        ['Other'],
        ['Multijet_MatrixMethod', proxies.MatrixMethod()]
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


class ChannelGenerator:

    def __init__(self, config, samples, decay_mode, signs, years, leptons, charges, btags, process_string, sample_config):
        self.config = config
        self.samples = samples
        self.decay_mode = decay_mode
        self.signs = signs
        self.years = years
        self.leptons = leptons
        self.charges = charges
        self.btags = btags
        self.process_string = process_string
        self.sample_config = sample_config

    def get_config(self):
        return self.config

    def make_channel(self, lumi, sign='', year='', lepton='', charge='', btag='', extra_rebin=1):
        channel_name = self.generate_channel_name(sign=sign, year=year, lepton=lepton, charge=charge, btag=btag)
        regions = self.generate_channel_regions(sign=sign, year=year, lepton=lepton, charge=charge, btag=btag)
        channel = {
            'regions': regions,
            'label': self.generate_channel_labels(sign=sign, lepton=lepton, charge=charge, btag=btag),
            'lumi': '+'.join(lumi),
            'extra_rebin': extra_rebin,
            'samples': [],
        }
        self.config['channels'][channel_name] = channel

        # print out
        print(f'Generated channel {channel_name} with {len(regions)} regions')
        for reg in regions:
            print(f'  {reg}')

        # add samples
        for sample in self.samples:
            if len(sample) == 1:
                channel['samples'] += [sample[0]]
            else:
                # generate proxy channel
                proxy = sample[1]
                proxy_channel_name = channel_name + "_" + proxy.get_name()
                proxy_channel = {
                    'make_plots': False,
                    'regions': proxy.format(regions),
                }
                self.config['channels'][proxy_channel_name] = proxy_channel
                channel['samples'] += [f'{sample[0]} | {proxy_channel_name}']

    def generate_channel_name(self, sign='', year='', lepton='', charge='', btag=''):
        if sign:
            sign = f'{sign}_'
        else:
            sign = 'OS-SS_'
        name = f'{sign}{format(year)}{format(lepton)}{format(charge)}{format(btag)}{self.decay_mode}'
        return name

    def generate_channel_regions(self, sign='', year='', lepton='', charge='', btag=''):
        regions = []
        if sign:
            regions = [f'{format(y)}{format(l)}{format(c)}{format(b)}{self.decay_mode}_{sign}'
                       for y in (self.years if not year else [year])
                       for l in (self.leptons if not lepton else [lepton])
                       for c in (self.charges if not charge else [charge])
                       for b in (self.btags if not btag else [btag])]
        else:
            regions = [f'{format(y)}{format(l)}{format(c)}{format(b)}{self.decay_mode}_OS'
                       for y in (self.years if not year else [year])
                       for l in (self.leptons if not lepton else [lepton])
                       for c in (self.charges if not charge else [charge])
                       for b in (self.btags if not btag else [btag])] + \
                [f'-{format(y)}{format(l)}{format(c)}{format(b)}{self.decay_mode}_SS'
                 for y in (self.years if not year else [year])
                 for l in (self.leptons if not lepton else [lepton])
                 for c in (self.charges if not charge else [charge])
                 for b in (self.btags if not btag else [btag])]
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
            row2 += f', {btag}'
        labels += [row2]
        labels += [label_dict[self.sample_config]]
        return labels
