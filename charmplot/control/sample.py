import logging
import ROOT
import sys

ROOT.gROOT.SetBatch(True)
logger = logging.getLogger(__name__)


class Sample(object):

    name = ""
    channel = ""
    add = []
    subtract = []
    color_scheme = None

    def __init__(self, name, color_scheme=None, **kwargs):
        # sample name
        self.name = name
        self.shortName = name

        # addition
        if 'add' in kwargs:
            add = kwargs.pop('add')
            if type(add) == list:
                self.add = add
            elif type(add) == str:
                self.add = [add]
        else:
            logger.critical("property 'add' missing in sample")
            sys.exit(1)

        # subtraction
        if 'subtract' in kwargs:
            subtract = kwargs.pop('subtract')
            if type(subtract) == list:
                self.subtract = subtract
            elif type(subtract) == str:
                self.subtract = [subtract]
        else:
            logger.critical("property 'subtract' missing in sample")
            sys.exit(1)

        # sample configuration
        if 'lineColor' in kwargs:
            self.lineColor = eval(kwargs.pop('lineColor'))
        else:
            self.lineColor = None
        if 'fillColor' in kwargs:
            self.fillColor = eval(kwargs.pop('fillColor'))
        else:
            self.fillColor = None

        if 'legendLabel' in kwargs:
            self.legendLabel = kwargs.pop('legendLabel')

        if 'shortName' in kwargs:
            self.shortName = kwargs.pop('shortName')

    def set_color_scheme(self, scheme):
        self.color_scheme = scheme

    def set_channel(self, channel):
        self.channel = channel

    def get_all(self):
        return self.add + self.subtract

    def __repr__(self):
        string = "+".join(self.add)
        for c in self.subtract:
            string += "-%s" % c
        return string
