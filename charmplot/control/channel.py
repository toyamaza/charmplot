

class Channel(object):

    name = ""
    label = ""
    lumi = ""
    add = []
    subtract = []

    def __init__(self, name, label, lumi, add, subtract=[]):
        self.name = name
        self.label = label
        self.lumi = lumi
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

    def get_all(self):
        return self.add + self.subtract

    def __repr__(self):
        string = "+".join(self.add)
        for c in self.subtract:
            string += "-%s" % c
        return string
