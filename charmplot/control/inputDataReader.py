import logging
import ROOT


from charmplot.common import utils


logger = logging.getLogger(__name__)


def skip_truth_channel(sample, c, input_file):
    if sample.ignoreTruthMatching:
        truth_channel = False
        if not (c.endswith("OS") or c.endswith("SS")):
            truth_channel = True
            logging.info(f"Truth channel detected: {c}")
        if truth_channel and input_file in sample.ignoreTruthMatching:
            logging.info(f"Skipping for file: {input_file}")
            return True
        if not truth_channel and input_file not in sample.ignoreTruthMatching:
            logging.info(f"Skipping for file: {input_file}")
            return True
    return False


class InputDataReader(object):

    # one input file for each sample
    input_files = {}
    channel_scale_factors = None

    def __init__(self, conf=None):
        self.config = conf
        self.read_input_files()

    def get_histogram_from_file(self, f, channel, sample, variable, c, extra_rebin=1, sys=None):
        logger.debug(f"In get_histogram_from_file with {f.GetName()} {channel.name} {sample.name}")
        # logger.debug(f"In get_histogram_from_file with {channel.name} {sample.name} {variable} {c} {f} extra_rebin: {extra_rebin}")
        h_name = "__".join([c, variable.name])
        h = None
        if sys:
            h_name_sys = "_-_".join([sys, "__".join([c, variable.name])])
            h = f.Get(h_name_sys)
            if not h:
                logger.warning(f"Systematic histogram with name {h_name_sys} not found! Using nominal instead.")
        if not h:
            h = f.Get(h_name)
        if not h:
            logger.warning(f"Histogram {h_name} not found for sample {sample.name} in file {f.GetName()}!")
            return None
        if sys:
            h = h.Clone(f"{h.GetName()}_{sys}_{f.GetName().replace('.root', '')}_{channel.name}")
        else:
            h = h.Clone(f"{h.GetName()}_{f.GetName().replace('.root', '')}_{channel.name}")
        h_new = utils.rebin_histogram(h, variable, extra_rebin)
        logger.info(f"Got histogram {h_new}")

        # scale histogram
        if self.channel_scale_factors and sample.shortName in self.channel_scale_factors:
            key = sample.shortName
            logging.info(f"Scaling histogram by {self.channel_scale_factors[key]}")
            h_new.Scale(self.channel_scale_factors[key])

        # scaling from MC config
        if sample.scaleMC and 'data' not in f.GetName():
            h_new.Scale(sample.scaleMC)
        return h_new

    def eff_divide(self, h_num, h_den):
        quot = h_num.Clone()
        quot.Divide(h_den)
        for i in range(h_num.GetNbinsX()):
            k = h_num.GetBinContent(i)
            n = h_den.GetBinContent(i)
            if n > 0 and n > k:
                print("error bits")
                print(str(k) + ", " + str(n))
                print(pow(k * (n - k) / (n**3), .5))
                quot.SetBinError(i, pow(k * (n - k) / (n**3), .5))
            else:
                quot.SetBinError(i, 0.0)
        return quot

    def get_histogram(self, sample, channel, variable, force_positive=False,
                      sys=None, suffix=None, integral_OS=None, integral_SS=None):
        logging.info(f":::get_histogram::: {sample.name} {channel.name}")
        h_plus = None
        h_minus = None
        h_total = None
        h_denom = None

        # scale factors for this channel
        self.channel_scale_factors = channel.scale_factors

        # extra rebin
        extra_rebin = channel.extra_rebin

        # take channels from the sample instead of global
        if sample.channel:
            channel = self.config.get_channel(sample.channel)
            logging.info(f"Using channel {channel.name} for sample {sample.name}")

        # looper over input files
        for input_file in sample.get_all():
            # check that fhe file exists
            if input_file not in self.input_files:
                logger.critical("No input file found for sample %s", sample.name)
                raise IOError("No input file found for sample %s", sample.name)
            f = self.input_files[input_file]
            for c in channel.add:
                logger.info(f"channel.add: {c} {sample.shortName}")
                if skip_truth_channel(sample, c, input_file):
                    continue
                h = self.get_histogram_from_file(f, channel, sample, variable, c, extra_rebin, sys)
                if not h:
                    continue
                if not h_plus:
                    h_plus = h.Clone("%s_%s_%s_plus" % (sample.name, channel.name, variable.name))
                else:
                    h_plus.Add(h)
            for c in channel.subtract:
                logger.info(f"channel.subtract: {c} {sample.shortName}")
                if skip_truth_channel(sample, c, input_file):
                    continue
                h = self.get_histogram_from_file(f, channel, sample, variable, c, extra_rebin, sys)
                if not h:
                    continue
                if not h_minus:
                    h_minus = h.Clone("%s_%s_%s_minus" % (sample.name, channel.name, variable.name))
                else:
                    h_minus.Add(h)
            for c in channel.divide:
                print("Division occurs")
                logger.info(f"channel.divide: {c} {sample.shortName}")
                if skip_truth_channel(sample, c, input_file):
                    continue
                h = self.get_histogram_from_file(f, channel, sample, variable, c, extra_rebin)
                if not h:
                    continue
                if not h_denom:
                    h_denom = h.Clone("%s_%s_%s" % (sample.name, channel.name, variable.name))
                else:
                    h_denom.Add(h)

        # need the plus histogram at this point...
        if not h_plus:
            return None

        # scale plus and minus if requested
        if integral_OS != None and h_plus:
            h_plus.Scale(integral_OS / h_plus.GetSumOfWeights())
        if integral_SS != None and h_minus:
            h_minus.Scale(integral_SS / h_minus.GetSumOfWeights())

        # add plus and minus
        h_total = h_plus.Clone("%s_%s_%s" % (sample.name, channel.name, variable.name))
        if h_minus:
            h_total.Add(h_minus, -1)

        if sample.offset:
            for i in range(0, h_total.GetNbinsX() + 2):
                h_total.SetBinContent(i, sample.offset)
        if sample.fit and h_total:
            h_total = utils.fit_histogram(h_total, sample.fit)
        if force_positive and h_total:
            # if 'MatrixMethod' not in sample.name:
            utils.set_to_positive(h_total, sys)
        if not sample.statError:
            logger.info(f"Set stat error of {sample} to zero.")
            for i in range(0, h_total.GetNbinsX() + 2):
                h_total.SetBinError(i, 0)
        if h_denom:
            h_total = self.eff_divide(h_total, h_denom)
        if suffix:
            h_total.SetName(f"{h_total.GetName()}_{suffix}")
            h_total.SetTitle(f"{h_total.GetName()}_{suffix}")
        logging.info(f":::END get_histogram::: {sample.name} {channel.name} {h_total}")
        return h_total

    def get_integral(self, sample, channel, variable, sys=None):
        assert len(channel.divide) == 0, "not implemented"
        integral_plus = 0
        integral_minus = 0
        if sample.channel:
            channel = self.config.get_channel(sample.channel)
            logging.info(f"Using channel {channel.name} for sample {sample.name}")
        for input_file in sample.get_all():
            if input_file not in self.input_files:
                logger.critical("No input file found for sample %s", sample.name)
                raise IOError("No input file found for sample %s", sample.name)
            f = self.input_files[input_file]
            for c in channel.add:
                logger.info(f"channel.add: {c} {sample.shortName}")
                if skip_truth_channel(sample, c, input_file):
                    continue
                h = self.get_histogram_from_file(f, channel, sample, variable, c, extra_rebin=1, sys=sys)
                if not h:
                    continue
                integral_plus += h.GetSumOfWeights()
            for c in channel.subtract:
                logger.info(f"channel.subtract: {c} {sample.shortName}")
                if skip_truth_channel(sample, c, input_file):
                    continue
                h = self.get_histogram_from_file(f, channel, sample, variable, c, extra_rebin=1, sys=sys)
                if not h:
                    continue
                integral_minus += h.GetSumOfWeights()
        return integral_plus, integral_minus

    def find_variables(self, sample, channel):
        all_vars = set()
        for input_file in sample.get_all():
            if input_file not in self.input_files:
                logger.critical("No input file found for sample %s: %s" % (sample.name, input_file))
                raise IOError("No input file found for sample %s: %s" % (sample.name, input_file))
            f = self.input_files[input_file]
            variables = {}
            for c in channel.get_all():
                logger.info(f"channel: {c}")
                for key in f.GetListOfKeys():
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
                if variables[var] > 0:
                    all_vars.add(var)
        return all_vars

    def read_input_file(self, sample):
        for input_file in sample.get_all():
            if input_file in self.input_files:
                continue
            f = ROOT.TFile("%s.root" % input_file, "READ")
            if not f:
                logger.critical("No input file found for sample %s", sample.name)
                raise IOError("No input file found for sample %s", sample.name)
            else:
                logger.info("Successfully read input file for sample %s", sample.name)
                self.input_files[input_file] = f

    def read_input_files(self):
        for s in self.config.get_data_and_mc():
            self.read_input_file(s)
