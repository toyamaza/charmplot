
from charmplot.control import channel
from charmplot.control import colorScheme
from charmplot.control import sample
from charmplot.control import tools
from charmplot.control import variable
import copy
import logging
import os


logger = logging.getLogger(__name__)


def parse_regions(regions):
    plus = []
    minus = []
    if type(regions) == list:
        for c in regions:
            if c.startswith("-"):
                minus += [c[1:]]
            else:
                plus += [c]
    elif type(regions) == str:
        plus += [regions]
    return plus, minus


class GlobalConfig(object):

    required_arguments = [
        'channels',
    ]

    # global samples configuration
    samples_config = {}

    data = None
    samples = []
    channels = []

    def get_data_and_mc(self):
        if self.data:
            return [self.data] + self.samples
        else:
            return self.samples

    def get_data(self):
        return self.data

    def get_bottom_mc(self):
        return self.default_samples[-1]

    def get_mc(self):
        return self.default_samples

    def get_sample_names(self):
        if hasattr(self, 'default_samples'):
            return [s.name for s in self.default_samples]
        else:
            return []

    def get_variables(self):
        return [self.get_var(v) for v in self.variables]

    def get_variable_names(self):
        return self.variables

    def get_var(self, v):
        if v not in self.variables:
            self.variables[v] = variable.Variable(v)
        return self.variables[v]

    def get_sample(self, name):
        for s in self.samples:
            if s.name == name:
                return s

    def get_channel(self, name):
        for c in self.channels:
            if c.name == name:
                return c

    def construct_sample(self, name):
        split_name = [x.strip() for x in name.split("|")]
        samp = self.samples_config[split_name[0]]
        s = sample.Sample(name, **samp)
        if len(split_name) == 2:
            s.set_channel(split_name[1])
        s.set_color_scheme(self.style)
        logger.debug(f"adding sample {s.name}: {s}")
        return s

    def read_variables(self, conf):
        variables = {}
        for name in self.variables_config:
            var = self.variables_config[name]
            v = variable.Variable(name, **var)
            variables[name] = v
        self.variables = variables

        # channel specific variables
        if 'channelsVariables' in conf:
            for name in conf['channelsVariables']:
                var = conf['channelsVariables'][name]
                v = variable.Variable(name, **var)
                self.variables[name] = v

    def read_samples(self, conf):
        samples = []
        if 'samples' in conf:
            for name in conf['samples']:
                s = self.construct_sample(name)
                samples += [s]
            self.samples = samples
            self.default_samples = copy.deepcopy(samples)

    def read_channel(self, conf):
        channels = []
        for name, val in conf['channels'].items():

            # parse inputs
            label = val['label'] if 'label' in val else ''
            lumi = val['lumi'] if 'lumi' in val else 0
            plus, minus = parse_regions(val['regions'])

            # make channel object
            chan = channel.Channel(name, label, lumi, plus, minus)
            channels += [chan]

            # append properties
            if 'samples' in val:
                chan.set_samples(val['samples'])
                for name in chan.samples:
                    if name not in self.get_sample_names():
                        s = self.construct_sample(name)
                        self.samples += [s]
            if 'make_plots' in val:
                chan.set_make_plots(val['make_plots'])
            if 'qcd_template' in val:
                chan.set_qcd_template(val['qcd_template'])
            if 'likelihood_fit' in val:
                chan.set_likelihood_fit(val['likelihood_fit'])
            if 'scale_factors' in val:
                chan.set_scale_factors(val['scale_factors'])
            if 'extra_rebin' in val:
                chan.set_extra_rebin(val['extra_rebin'])
            if 'mass_fit' in val:
                chan.set_mass_fit(val['mass_fit'])
        self.channels = channels

    def parse_confing(self, conf):
        for arg in self.required_arguments:
            if arg not in conf:
                logger.critical("argument %s not found in config!" % arg)
                raise Exception("Invalid config")

        # read data
        if 'data' in conf:
            samp = self.samples_config[conf['data']]
            self.data = sample.Sample(conf['data'], **samp)

        # read variables
        self.read_variables(conf)

        # read samples
        self.read_samples(conf)

        # read channels
        self.read_channel(conf)

    def __init__(self, config_name):
        self.config_name = config_name

        # analysis specific config
        conf = tools.parse_yaml(config_name)

        # global samples config
        assert 'samplesConf' in conf
        self.samples_config = tools.parse_yaml(os.path.join('samples', conf['samplesConf']))

        # color scheme
        color_scheme = getattr(colorScheme, self.samples_config['colorScheme'])
        self.style = color_scheme()

        # global variables config
        assert 'variablesConf' in conf
        self.variables_config = tools.parse_yaml(os.path.join('variables', conf['variablesConf']))

        # parse config
        self.parse_confing(conf)
