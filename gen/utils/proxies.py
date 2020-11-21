#!/usr/bin/env python
import abc


class ProxyChannel:
    @property
    @abc.abstractmethod
    def name(self):
        pass

    @abc.abstractmethod
    def format(self, regions):
        pass

    def get_name(self):
        return self.name


class MatrixMethod(ProxyChannel):
    name = 'MatrixMethod'

    def format(self, regions):
        out = []
        for reg in regions:
            sign = ""
            anti_sign = "-"
            if reg.startswith("-"):
                reg = reg[1:]
                sign = "-"
                anti_sign = ""
            out += [f"{sign}AntiTight_{reg}"]
            out += [f"{anti_sign}Tight_{reg}"]
        return out


class Matched(ProxyChannel):
    name = 'Matched'

    def format(self, regions):
        return [reg + "_Matched" for reg in regions]


class NoMatch(ProxyChannel):
    name = 'NoMatch'

    def format(self, regions):
        return [reg + "_Other" for reg in regions]


class Rest(ProxyChannel):
    name = 'Rest'

    def format(self, regions):
        out = []
        for reg in regions:
            anti_sign = "-"
            if reg.startswith("-"):
                reg = reg[1:]
                anti_sign = ""
            out += [f"{anti_sign}{reg}_Matched"]
            out += [f"{anti_sign}{reg}_Other"]
        return regions + out
