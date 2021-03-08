#!/usr/bin/env python
import abc


class ProxyChannel:
    os_minus_ss_fit_configuration = False

    def __init__(self, os_minus_ss_fit_configuration: bool = False, loose_sr: bool = False):
        self.os_minus_ss_fit_configuration = os_minus_ss_fit_configuration
        self.loose_sr = loose_sr

    def format(self, regions):
        out = []
        if not self.os_minus_ss_fit_configuration:
            out = regions
        else:
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
        if not self.loose_sr:
            return out
        else:
            antiSR = []
            for reg in out:
                antiSR += [reg.replace("_SR", "_Anti_SR")]
            return antiSR + out

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

    def __init__(self, os_minus_ss_fit_configuration: bool = False, loose_sr: bool = False, subtract_mj: bool = False):
        super().__init__(os_minus_ss_fit_configuration=os_minus_ss_fit_configuration, loose_sr=loose_sr)
        self.subtract_mj = subtract_mj

    def get_regions(self, regions):
        out = []
        for reg in regions:
            sign = ""
            anti_sign = "-"
            if reg.startswith("-"):
                reg = reg[1:]
                sign = "-"
                anti_sign = ""
            if "_OS" in reg:
                reg_SS = reg.replace("_OS", "_SS")
                out += [reg_SS]
                if self.subtract_mj:
                    out += [f"{anti_sign}AntiTight_{reg_SS}"]
                    out += [f"{sign}Tight_{reg_SS}"]
            else:
                out += [reg]
                if self.subtract_mj:
                    out += [f"{anti_sign}AntiTight_{reg}"]
                    out += [f"{sign}Tight_{reg}"]
        return self.format(out)


class MatrixMethod(ProxyChannel):
    name = 'MatrixMethod'

    def __init__(self, os_minus_ss_fit_configuration: bool = False, loose_sr: bool = False, fake_factor: bool = False):
        super().__init__(os_minus_ss_fit_configuration=os_minus_ss_fit_configuration, loose_sr=loose_sr)
        self.fake_factor = fake_factor
        if self.fake_factor:
            self.name = 'FakeFactor'

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
            if not self.fake_factor:
                out += [f"{anti_sign}Tight_{reg}"]
        return self.format(out)


class Matched(ProxyChannel):
    name = 'Matched'

    def get_regions(self, regions):
        return self.format([reg + "_Matched" for reg in regions])


class GenericChannel(ProxyChannel):
    name = ""

    def __init__(self, os_minus_ss_fit_configuration: bool = False, loose_sr: bool = False, region: str = "", name: str = "", regions_override: list = []):
        super().__init__(os_minus_ss_fit_configuration=os_minus_ss_fit_configuration, loose_sr=loose_sr)
        self.region = region
        self.regions_override = regions_override
        if name == "":
            self.name = region
        else:
            self.name = name

    def get_regions(self, regions):
        if self.regions_override:
            return self.format(self.regions_override)
        if type(self.region) == str:
            return self.format([f"{reg}_{self.region}" for reg in regions])
        else:
            return self.format([f"{reg1}_{reg2}" for reg1 in regions for reg2 in self.region])


class MisMatched(ProxyChannel):
    name = 'MisMatched'
    truth_name = 'MisMatched'

    def __init__(self, os_minus_ss_fit_configuration: bool = False, loose_sr: bool = False, pdgId: str = ""):
        super().__init__(os_minus_ss_fit_configuration=os_minus_ss_fit_configuration, loose_sr=loose_sr)
        self.pdgId = pdgId
        if self.pdgId != "":
            self.name = f"{self.pdgId}{self.name}"
            self.truth_name = self.name
        if loose_sr:
            self.name += "_Loose"

    def get_regions(self, regions):
        return self.format([reg + f"_{self.truth_name}" for reg in regions])


class MatchedFid(ProxyChannel):
    name = 'MatchedFid'

    def get_regions(self, regions):
        return self.format([reg + "_MatchedFid" for reg in regions])


class MatchedNoFid(ProxyChannel):
    name = 'MatchedNoFid'

    def get_regions(self, regions):
        return self.format([reg + "_MatchedNoFid" for reg in regions])


class NoMatch(ProxyChannel):
    name = 'NoMatch'

    def get_regions(self, regions):
        return self.format([reg + "_Other" for reg in regions])


class NoMatchBackground(ProxyChannel):
    name = 'NoMatchBackground'

    def __init__(self, os_minus_ss_fit_configuration: bool = False, loose_sr: bool = False):
        super().__init__(os_minus_ss_fit_configuration=os_minus_ss_fit_configuration, loose_sr=loose_sr)
        if loose_sr:
            self.name += "_Loose"

    def get_regions(self, regions):
        return self.format([reg + "_Other" for reg in regions] +
                           [reg + "_HardMisMatched" for reg in regions])


class MatchedCharm(ProxyChannel):
    name = 'MatchedCharm'

    def __init__(self, os_minus_ss_fit_configuration: bool = False, loose_sr: bool = False):
        super().__init__(os_minus_ss_fit_configuration=os_minus_ss_fit_configuration, loose_sr=loose_sr)
        if loose_sr:
            self.name += "_Loose"

    def get_regions(self, regions):
        return self.format([reg + "_431MisMatched" for reg in regions] +
                           [reg + "_413MisMatched" for reg in regions] +
                           [reg + "_421MisMatched" for reg in regions] +
                           [reg + "_BaryonMisMatched" for reg in regions])
# [reg + "_MisMatched" for reg in regions] +
#                            [reg + "_MatchedNoFid" for reg in regions] +


class MatchedDplus(ProxyChannel):
    name = 'MatchedDplus'

    def get_regions(self, regions):
        return self.format([reg + "_Matched" for reg in regions] +
                           [reg + "_411MisMatched" for reg in regions])


class Rest(ProxyChannel):
    name = 'Rest'

    def __init__(self, os_minus_ss_fit_configuration: bool = False, loose_sr: bool = False, allowed_regions: list = [],
                 name: str = "", exclude_mismatched: bool = False, include_other: bool = False):
        super().__init__(os_minus_ss_fit_configuration=os_minus_ss_fit_configuration, loose_sr=loose_sr)
        self.allowed_regions = allowed_regions
        if name:
            if not name.startswith("_"):
                name = "_" + name
            self.name += name
        self.exclude_mismatched = exclude_mismatched
        self.include_other = include_other
        if loose_sr:
            self.name += "_Loose"

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
            if not self.include_other:
                out += [f"{anti_sign}{reg}_Other"]
            if self.exclude_mismatched:
                out += [f"{anti_sign}{reg}_411MisMatched"]
                out += [f"{anti_sign}{reg}_413MisMatched"]
                out += [f"{anti_sign}{reg}_421MisMatched"]
                out += [f"{anti_sign}{reg}_431MisMatched"]
                out += [f"{anti_sign}{reg}_BaryonMisMatched"]
                out += [f"{anti_sign}{reg}_MatchedNoFid"]
                out += [f"{anti_sign}{reg}_MisMatched"]
        return self.format(regs + out)
