import logging


logger = logging.getLogger(__name__)


class Variable(object):

    label = None
    unit = None
    x_range = None
    ratio_range = None
    stage_out = True
    make_plots = True
    allow_rebin = True
    rebin = 1
    do_overflow = True
    do_underflow = True
    dstar_tail_rebin = False
    bins = None
    logx = False
    xbins = None
    name_override = None

    def __init__(self, name, **kwargs):
        self.name = name
        for key in kwargs:
            logger.debug(f"setting {key} to {kwargs[key]}")
            setattr(self, key, kwargs[key])

        if self.name_override:
            self.name = self.name_override
