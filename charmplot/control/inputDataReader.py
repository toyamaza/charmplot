import logging
import os
import ROOT


from charmplot.common import utils


logger = logging.getLogger(__name__)


class InputDataReader(object):

    # one input file for each sample
    input_files = {}
    channel_scale_factors = None

    def __init__(self, conf=None):
        self.config = conf
        self.read_input_files()

    def get_histogram_from_file(self, f, channel, sample, variable, c):
        h_name = os.path.join(c, "__".join([c, variable.name]))
        h = f.Get(h_name)
        if not h:
            logger.warning("Histogram %s for variable %s not found in channel %s for sample %s in file %s" % (
                h_name, variable.name, c, sample.name, f))
            return None
        h = h.Clone(h.GetName() + "_temp")
        utils.rebin_histogram(h, variable, channel.extra_rebin)
        utils.set_to_positive(h)

        # scale histogram
        scale_factors = utils.read_scale_factors(self.channel_scale_factors)
        self.scale_histogram(h, sample, c, scale_factors)

        return h

    def scale_histogram(self, h, sample, channel, scale_factors):
        # scale histogram if given input scale factors
        if scale_factors:
            sf = [1., 0.]
            logger.info(f"s.shortName: {sample.shortName} {channel}")
            for s in self.channel_scale_factors['scale_factors']:
                sample_channel = s.split(" | ")
                if len(sample_channel) == 2 and sample_channel[1] not in channel:
                    continue
                if sample.shortName == sample_channel[0]:
                    if self.channel_scale_factors['scale_factors'][s] in scale_factors.keys():
                        sf = scale_factors[self.channel_scale_factors['scale_factors'][s]]
            h.Scale(sf[0])
            logger.info(f"Scaling histogram {sample.name} in {channel} by {sf[0]}")

    def get_histogram(self, sample, channel, variable):
        h_total = None
        # scale factors for this channel
        self.channel_scale_factors = channel.scale_factors
        if sample.channel:
            logging.info(f"Using channel {sample.channel} for sample {sample.name}")
            channel = self.config.get_channel(sample.channel)
        for input_file in sample.get_all():
            if input_file in sample.add:
                weight = 1.
            elif input_file in sample.subtract:
                weight = -1.
            if input_file not in self.input_files:
                logger.critical("No input file found for sample %s", sample.name)
                raise IOError("No input file found for sample %s", sample.name)
            f = self.input_files[input_file]
            for c in channel.add:
                logger.info(f"channel.add: {c} {sample.shortName}")
                h = self.get_histogram_from_file(f, channel, sample, variable, c)
                if not h:
                    continue
                if not h_total:
                    h_total = h.Clone("%s_%s_%s" % (sample.name, channel.name, variable.name))
                else:
                    h_total.Add(h, weight)
            for c in channel.subtract:
                logger.info(f"channel.subtract: {c} {sample.shortName}")
                h = self.get_histogram_from_file(f, channel, sample, variable, c)
                if not h:
                    continue
                if not h_total:
                    return None
                h_total.Add(h, -1 * weight)
        return h_total

    def find_variables(self, sample, channel):
        all_vars = set()
        for input_file in sample.get_all():
            if input_file not in self.input_files:
                logger.critical("No input file found for sample %s: %s" % (sample.name, input_file))
                raise IOError("No input file found for sample %s: %s" % (sample.name, input_file))
            f = self.input_files[input_file]
            variables = {}
            for c in channel.get_all():
                directory = f.Get(c)
                for key in directory.GetListOfKeys():
                    name = key.GetName()
                    var_name = name
                    if name.startswith(c):
                        var_name = var_name[len(c):]
                        while var_name.startswith("_"):
                            var_name = var_name[1:]
                    if var_name not in variables:
                        variables[var_name] = 1
                    else:
                        variables[var_name] += 1
            for var in variables:
                if variables[var] == len(channel.get_all()):
                    all_vars.add(var)
        return all_vars

    def read_input_files(self):
        for s in self.config.get_data_and_mc():
            for input_file in s.get_all():
                if input_file in self.input_files:
                    continue
                f = ROOT.TFile("%s.root" % input_file, "READ")
                if not f:
                    logger.critical("No input file found for sample %s", s.name)
                    raise IOError("No input file found for sample %s", s.name)
                else:
                    logger.info("Successfully read input file for sample %s", s.name)
                    self.input_files[input_file] = f
