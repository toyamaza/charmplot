
from charmplot.control import channel
from charmplot.control import colorScheme
from charmplot.control import sample
from charmplot.control import tools
from charmplot.control import variable
import logging
import os


logger = logging.getLogger(__name__)


class GlobalConfig(object):

    required_arguments = [
        'channels',
    ]

    # global samples configuration
    samples_config = {}

    data = None
    samples = []
    channels = []

    def get_bottom_mc(self):
        return self.samples[-1]

    def get_mc(self):
        return self.samples

    def get_sample_names(self):
        return [s.name for s in self.samples]

    def get_data(self):
        return self.data

    def get_data_and_mc(self):
        if self.data:
            return [self.data] + self.samples
        else:
            return self.samples

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
        samp = self.samples_config[name]
        s = sample.Sample(name, **samp)
        s.set_color_scheme(self.style)
        logger.debug(f"adding sample {s.name}: {s}")
        return s

    def parse_confing(self, conf):
        for arg in self.required_arguments:
            if arg not in conf:
                logger.critical("argument %s not found in config!" % arg)
                raise Exception("Invalid config")

        # color scheme
        color_scheme = getattr(colorScheme, conf['colorScheme'])
        self.style = color_scheme()

        # read data
        if 'data' in conf:
            samp = self.samples_config[conf['data']]
            self.data = sample.Sample(conf['data'], **samp)

        # read variables
        variables = {}
        for name in self.variables_config:
            var = self.variables_config[name]
            v = variable.Variable(name, **var)
            variables[name] = v
        self.variables = variables

        # read samples
        samples = []
        if 'samples' in conf:
            for name in conf['samples']:
                s = self.construct_sample(name)
                samples += [s]
            self.samples = samples

        # read channels
        channels = []
        for name, val in conf['channels'].items():
            regions = val['regions']
            label = val['label'] if 'label' in val else ''
            lumi = val['lumi'] if 'lumi' in val else 0
            if type(lumi) == int:
                lumi = str(lumi)
            if type(regions) == list:
                plus = []
                minus = []
                for c in regions:
                    if c.startswith("-"):
                        minus += [c[1:]]
                    else:
                        plus += [c]
                chan = channel.Channel(name, label, lumi, plus, minus)
            elif type(regions) == str:
                chan = channel.Channel(name, label, lumi, regions)
            else:
                logger.critical("unrecognized channel type for " % name)
                raise Exception("Invalid config")
            if 'samples' in val:
                chan.set_samples(val['samples'])
                for name in chan.samples:
                    if not name in self.get_sample_names():
                        s = self.construct_sample(name)
                        self.samples += [s]
            if 'qcd_template' in val:
                chan.set_qcd_template(val['qcd_template'])
            channels += [chan]
        self.channels = channels

    def __init__(self, config_name):
        # analysis specific config
        conf = tools.parse_yaml(config_name)

        # global samples config
        assert 'samplesConf' in conf
        self.samples_config = tools.parse_yaml(os.path.join('samples', conf['samplesConf']))

        # global variables config
        assert 'variablesConf' in conf
        self.variables_config = tools.parse_yaml(os.path.join('variables', conf['variablesConf']))

        # parse config
        self.parse_confing(conf)
