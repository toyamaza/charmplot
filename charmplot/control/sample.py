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
    offset = None
    color_scheme = None
    statError = True
    scaleMC = None
    ghost = False
    fit = ""
    ignoreTruthMatching = []
    allowRebin = True
    ignoreForcePositive = False
    fitRebin = False

    def __init__(self, name, color_scheme=None, **kwargs):
        # sample name
        self.name = name
        self.shortName = name

        # offset
        if 'offset' in kwargs:
            self.offset = float(kwargs.pop('offset'))

        # addition
        if 'add' in kwargs:
            add = kwargs.pop('add')
            if type(add) == list:
                self.add = add
            elif type(add) == str:
                self.add = [add]
        elif not self.offset:
            logger.critical("property 'add' or 'offset' missing in sample")
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

        if 'lineStyle' in kwargs:
            self.lineStyle = int(kwargs.pop('lineStyle'))
        else:
            self.lineStyle = 1

        if 'fillColor' in kwargs:
            self.fillColor = eval(kwargs.pop('fillColor'))
        else:
            self.fillColor = None

        if 'legendLabel' in kwargs:
            self.legendLabel = kwargs.pop('legendLabel')

        if 'shortName' in kwargs:
            self.shortName = kwargs.pop('shortName')

        if 'statError' in kwargs:
            self.statError = kwargs.pop('statError')

        if 'scaleMC' in kwargs:
            if type(kwargs['scaleMC']) == str:
                self.scaleMC = float(eval(kwargs.pop('scaleMC')))
            else:
                self.scaleMC = float(kwargs.pop('scaleMC'))

        if 'fit' in kwargs:
            self.fit = str(kwargs.pop('fit'))

        if 'ignoreTruthMatching' in kwargs:
            self.ignoreTruthMatching = str(kwargs.pop('ignoreTruthMatching'))

        if 'ghost' in kwargs:
            self.ghost = bool(kwargs.pop('ghost'))

        if 'allowRebin' in kwargs:
            self.allowRebin = bool(kwargs.pop('allowRebin'))

        if 'ignoreForcePositive' in kwargs:
            self.ignoreForcePositive = bool(kwargs.pop('ignoreForcePositive'))

        if 'fitRebin' in kwargs:
            self.fitRebin = bool(kwargs.pop('fitRebin'))

    def set_color_scheme(self, scheme):
        self.color_scheme = scheme

    def set_channel(self, channel):
        self.channel = channel

    def get_all(self):
        return self.add + self.subtract

    def __repr__(self):
        return self.name
        # string = "+".join(self.add)
        # for c in self.subtract:
        #     string += "-%s" % c
        # return string
