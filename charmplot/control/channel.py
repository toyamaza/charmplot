

class Channel(object):

    qcd_template = None
    likelihood_fit = None
    mass_fit = None
    scale_factors = None
    make_plots = True
    save_to_file = True
    print_scale_factors = True
    force_positive = False
    extra_rebin = 1
    extra_scale = 1.0
    name = ""
    lumi = ""
    label = []
    add = []
    subtract = []
    divide = []
    samples = []
    trex_subtraction = {}
    replacement_samples = {}

    def __init__(self, name, label, lumi, add, subtract=[], divide = [], samples=[]):
        self.name = name
        self.lumi = lumi
        # label
        if type(label) == list:
            self.label = label
        elif type(label) == str:
            self.label = [label]
        # addition
        if type(add) == list:
            self.add = add
        elif type(add) == str:
            self.add = [add]
        # subtraction
        if type(subtract) == list:
            self.subtract = subtract
        elif type(subtract) == str:
            self.subtract = [subtract]
        # division
        if type(divide) == list:
            self.divide = divide
        elif type(divide) == str:
            self.divide = [divide]
        # samples
        if type(samples) == list:
            self.samples = samples
        elif type(samples) == str:
            self.samples = [samples]

    def set_make_plots(self, make_plots):
        self.make_plots = make_plots

    def set_qcd_template(self, qcd_template):
        self.qcd_template = qcd_template

    def set_samples(self, samples):
        self.samples = samples

    def set_likelihood_fit(self, fit):
        self.likelihood_fit = fit

    def set_scale_factors(self, scale_factors):
        self.scale_factors = scale_factors

    def set_extra_rebin(self, extra_rebin):
        self.extra_rebin = extra_rebin

    def set_extra_scale(self, extra_scale):
        self.extra_scale = extra_scale

    def set_mass_fit(self, mass_fit):
        self.mass_fit = mass_fit

    def set_save_to_file(self, save_to_file):
        self.save_to_file = save_to_file

    def set_print_scale_factors(self, print_scale_factors):
        self.print_scale_factors = print_scale_factors

    def set_force_positive(self, force_positive):
        self.force_positive = force_positive

    def set_trex_subtraction(self, trex_subtraction):
        self.trex_subtraction = trex_subtraction

    def set_replacement_samples(self, replacement_samples):
        self.replacement_samples = replacement_samples

    def get_all(self):
        return self.add + self.subtract

    def __repr__(self):
        string = "+".join(self.add)
        for c in self.subtract:
            string += "-%s" % c
        return string
