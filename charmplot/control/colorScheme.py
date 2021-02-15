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
        ROOT.TColor(10000, *RGB("#c1d2fc"), "10000"),  # Sapphire
        ROOT.TColor(10001, *RGB("#93aef7"), "10001"),
        ROOT.TColor(10002, *RGB("#698be9"), "10002"),
        ROOT.TColor(10003, *RGB("#4168d4"), "10003"),
        ROOT.TColor(10004, *RGB("#2049ae"), "10004"),
        ROOT.TColor(10005, *RGB("#002b89"), "10005"),
        ROOT.TColor(11000, *RGB("#89e5a5"), "11000"),  # Emerald
        ROOT.TColor(11001, *RGB("#48c576"), "11001"),
        ROOT.TColor(11002, *RGB("#1da256"), "11002"),
        ROOT.TColor(11003, *RGB("#007f37"), "11003"),
        ROOT.TColor(11004, *RGB("#005c1b"), "11004"),
        ROOT.TColor(11005, *RGB("#003c00"), "11005"),
        ROOT.TColor(12000, *RGB("#fccb38"), "12000"),  # Amber
        ROOT.TColor(12001, *RGB("#dfa400"), "12001"),
        ROOT.TColor(12002, *RGB("#b88300"), "12002"),
        ROOT.TColor(12003, *RGB("#936300"), "12003"),
        ROOT.TColor(12004, *RGB("#704400"), "12004"),
        ROOT.TColor(12005, *RGB("#502600"), "12005"),
        ROOT.TColor(13000, *RGB("#d1d1c5"), "13000"),  # Grey pearl
        ROOT.TColor(13001, *RGB("#aeaea0"), "13001"),
        ROOT.TColor(13002, *RGB("#8c8c7f"), "13002"),
        ROOT.TColor(13003, *RGB("#6b6b5e"), "13003"),
        ROOT.TColor(13004, *RGB("#4c4c40"), "13004"),
        ROOT.TColor(13005, *RGB("#2f2f24"), "13005"),
        ROOT.TColor(14000, *RGB("#fac6bf"), "14000"),  # Ruby
        ROOT.TColor(14001, *RGB("#fa9488"), "14001"),
        ROOT.TColor(14002, *RGB("#e36a62"), "14002"),
        ROOT.TColor(14003, *RGB("#c6403d"), "14003"),
        ROOT.TColor(14004, *RGB("#9d211f"), "14004"),
        ROOT.TColor(14005, *RGB("#720000"), "14005"),
        ROOT.TColor(15000, *RGB("#ecc2fd"), "15000"),  # Amethyst
        ROOT.TColor(15001, *RGB("#cf98ed"), "15001"),
        ROOT.TColor(15002, *RGB("#a975d4"), "15002"),
        ROOT.TColor(15003, *RGB("#8453b7"), "15003"),
        ROOT.TColor(15004, *RGB("#623394"), "15004"),
        ROOT.TColor(15005, *RGB("#3f1472"), "15005"),
        ROOT.TColor(16000, *RGB("#dfff4f"), "16000"),  # Odd colors
        ROOT.TColor(17000, *RGB("#68FF00"), "17000"),
        ROOT.TColor(18000, *RGB("#13F4EF"), "18000"),
        ROOT.TColor(19000, *RGB("#FF005C"), "19000"),
        ROOT.TColor(20000, *RGB("#006FFF"), "20000"),
        ROOT.TColor(21000, *RGB("#232b2b"), "21000"),
    ]
    return scheme
