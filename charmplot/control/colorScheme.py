import ROOT
import logging

logger = logging.getLogger(__name__)

# Colors taken from:
# http://color.adobe.com/create/color-wheel
# http://colorbrewer2.org/
# https://www.vis4.net/blog/2013/09/mastering-multi-hued-color-scales/


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


def cheer_up_emo_kid_extended():
    logger.info("Setting color scheme to 'cheer_up_emo_kid'")
    scheme = [ROOT.TColor(10000 + i, *RGB(c)) for i, c in enumerate([
        "#8fffff",
        "#4ecdc4",
        "#00948c",
        "#005d58",
        "#FF6B6B",
        "#C44D58",
        "#C7F464",
        "#C9C9C9",
        "#AE00FF",
    ])]
    return scheme


def gemstone_hues():
    logger.info("Setting color scheme to 'gemstone_hues'")
    scheme = [
        ROOT.TColor(10000, *RGB("#c1d2fc")),  # Sapphire
        ROOT.TColor(10001, *RGB("#93aef7")),
        ROOT.TColor(10002, *RGB("#698be9")),
        ROOT.TColor(10003, *RGB("#4168d4")),
        ROOT.TColor(10004, *RGB("#2049ae")),
        ROOT.TColor(10005, *RGB("#002b89")),
        ROOT.TColor(11000, *RGB("#89e5a5")),  # Emerald
        ROOT.TColor(11001, *RGB("#48c576")),
        ROOT.TColor(11002, *RGB("#1da256")),
        ROOT.TColor(11003, *RGB("#007f37")),
        ROOT.TColor(11004, *RGB("#005c1b")),
        ROOT.TColor(11005, *RGB("#003c00")),
        ROOT.TColor(12000, *RGB("#fccb38")),  # Amber
        ROOT.TColor(12001, *RGB("#dfa400")),
        ROOT.TColor(12002, *RGB("#b88300")),
        ROOT.TColor(12003, *RGB("#936300")),
        ROOT.TColor(12004, *RGB("#704400")),
        ROOT.TColor(12005, *RGB("#502600")),
        ROOT.TColor(13000, *RGB("#d1d1c5")),  # Grey pearl
        ROOT.TColor(13001, *RGB("#aeaea0")),
        ROOT.TColor(13002, *RGB("#8c8c7f")),
        ROOT.TColor(13003, *RGB("#6b6b5e")),
        ROOT.TColor(13004, *RGB("#4c4c40")),
        ROOT.TColor(13005, *RGB("#2f2f24")),
        ROOT.TColor(14000, *RGB("#fac6bf")),  # Ruby
        ROOT.TColor(14001, *RGB("#fa9488")),
        ROOT.TColor(14002, *RGB("#e36a62")),
        ROOT.TColor(14003, *RGB("#c6403d")),
        ROOT.TColor(14004, *RGB("#9d211f")),
        ROOT.TColor(14005, *RGB("#720000")),
        ROOT.TColor(15000, *RGB("#ecc2fd")),  # Amethyst
        ROOT.TColor(15001, *RGB("#cf98ed")),
        ROOT.TColor(15002, *RGB("#a975d4")),
        ROOT.TColor(15003, *RGB("#8453b7")),
        ROOT.TColor(15004, *RGB("#623394")),
        ROOT.TColor(15005, *RGB("#3f1472")),
    ]
    return scheme
