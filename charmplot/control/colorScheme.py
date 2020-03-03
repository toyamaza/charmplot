import ROOT
import logging

logger = logging.getLogger(__name__)

# Colors taken from:
# http://color.adobe.com/create/color-wheel
# http://colorbrewer2.org/


def RGB(string):
    return tuple(int(string[i + 1:i + 3], 16) / 255. for i in (0, 2, 4))


# https://www.colourlovers.com/palette/1930/cheer_up_emo_kid
def cheer_up_emo_kid():
    logger.info("Setting color scheme to 'cheer_up_emo_kid'")
    scheme = [ROOT.TColor(10000 + i, *RGB(c)) for i, c in enumerate([
        "#4ECDC4",
        "#556270",
        "#FF6B6B",
        "#C44D58",
        "#C7F464",
        "#C9C9C9",
        "#AE00FF",
    ])]
    return scheme
