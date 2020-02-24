import ROOT
import logging

logger = logging.getLogger(__name__)

# Colors taken from:
# http://color.adobe.com/create/color-wheel
# http://colorbrewer2.org/


def RGB(string):
    return tuple(int(string[i + 1:i + 3], 16) / 255. for i in (0, 2, 4))


def berkeley_wjets():
    logger.info("Setting color scheme to 'berkeley_wjets'")
    scheme = [ROOT.TColor(10000 + i, *RGB(c)) for i, c in enumerate([
        "#93c4d2",
        "#00429d",
        "#ffa59e",
        "#93003a",
        "#ffffbf",
        "#c9c9c9",
    ])]
    return scheme
