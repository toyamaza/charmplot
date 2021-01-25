#!/usr/bin/env python
import abc


class ProxyChannel:
    os_minus_ss_fit_configuration = False

    def __init__(self, os_minus_ss_fit_configuration: bool = False):
        self.os_minus_ss_fit_configuration = os_minus_ss_fit_configuration

    def format(self, regions):
        if not self.os_minus_ss_fit_configuration:
            return regions
        else:
            out = []
            for reg in regions:
                out += [reg]
                anti_sign = "-"
                reg_nosign = reg
                if reg.startswith("-"):
                    anti_sign = ""
                    reg_nosign = reg[1:]
                if "_OS" in reg_nosign:
                    reg_nosign = reg_nosign.replace("_OS", "_SS")
                else:
                    reg_nosign = reg_nosign.replace("_SS", "_OS")
                out += [f"{anti_sign}{reg_nosign}"]
            return out

    @property
    @abc.abstractmethod
    def name(self):
        pass

    @abc.abstractmethod
    def get_regions(self, regions):
        pass

    def get_name(self):
        return self.name


class PlainChannel(ProxyChannel):
    name = 'Plain'

    def get_regions(self, regions):
        return self.format(regions)


class MockMC(ProxyChannel):
    name = 'MockMC'

    def get_regions(self, regions):
        out = []
        for reg in regions:
            if "_OS" in reg:
                out += [reg.replace("_OS", "_SS")]
            else:
                out += [reg]
        return self.format(out)


class MatrixMethod(ProxyChannel):
    name = 'MatrixMethod'

    def get_regions(self, regions):
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
        return self.format(out)


class Matched(ProxyChannel):
    name = 'Matched'

    def get_regions(self, regions):
        return self.format([reg + "_Matched" for reg in regions])


class MatchedFid(ProxyChannel):
    name = 'MatchedFid'

    def get_regions(self, regions):
        return self.format([reg + "_MatchedFid" for reg in regions])


class MatchedNoFid(ProxyChannel):
    name = 'MatchedNoFid'

    def get_regions(self, regions):
        return self.format([reg + "_Matched" for reg in regions])


class NoMatch(ProxyChannel):
    name = 'NoMatch'

    def get_regions(self, regions):
        return self.format([reg + "_Other" for reg in regions])


# class Rest(ProxyChannel):
#     name = 'Rest'

#     def get_regions(self, regions):
#         out = []
#         for reg in regions:
#             anti_sign = "-"
#             if reg.startswith("-"):
#                 reg = reg[1:]
#                 anti_sign = ""
#             out += [f"{anti_sign}{reg}_Matched"]
#             out += [f"{anti_sign}{reg}_MatchedFid"]
#             out += [f"{anti_sign}{reg}_Other"]
#         return self.format(regions + out)

class Rest(ProxyChannel):
    name = 'Rest'

    def __init__(self, os_minus_ss_fit_configuration: bool = False, allowed_regions: list = [], name: str = ""):
        super().__init__(os_minus_ss_fit_configuration)
        self.allowed_regions = allowed_regions
        if name:
            if not name.startswith("_"):
                name = "_" + name
            self.name += name

    def get_regions(self, regions):
        out = []
        regs = []
        for reg in regions:
            pass_region = False
            for allowed_reg in self.allowed_regions:
                if allowed_reg in reg:
                    pass_region = True
                    break
            if len(self.allowed_regions) and not pass_region:
                continue
            regs += [reg]
            anti_sign = "-"
            if reg.startswith("-"):
                reg = reg[1:]
                anti_sign = ""
            out += [f"{anti_sign}{reg}_Matched"]
            out += [f"{anti_sign}{reg}_Other"]
        return self.format(regs + out)
