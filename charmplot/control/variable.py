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

    def __init__(self, name, **kwargs):
        self.name = name
        for key in kwargs:
            logger.debug(f"setting {key} to {kwargs[key]}")
            setattr(self, key, kwargs[key])
