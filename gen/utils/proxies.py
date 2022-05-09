#!/usr/bin/env python
import abc


class ProxyChannel:
    os_ss_sub = False

    def __init__(self, os_ss_sub: bool = False, loose_sr: bool = False, loose_sr_only: bool = False,
                 ss_only: bool = False, os_only: bool = False, force_1tag: bool = False, add_1tag: bool = False):
        self.os_ss_sub = os_ss_sub
        self.loose_sr = loose_sr
        self.loose_sr_only = loose_sr_only
        self.ss_only = ss_only
        self.os_only = os_only
        self.force_1tag = force_1tag
        self.add_1tag = add_1tag

        # inconsistent config
        assert not (ss_only and os_only)

    def format(self, regions):

        # change all 0-tag to 1-tag
        if self.force_1tag:
            for i in range(len(regions)):
                if "0tag" in regions[i]:
                    regions[i] = regions[i].replace("0tag", "1tag")

        # add 1-tag in addition to 0-tag
        elif self.add_1tag:
            additional = []
            for reg in regions:
                if "0tag" in reg:
                    additional += [reg.replace("0tag", "1tag")]
            regions += additional

        out = []
        if not self.os_ss_sub:
            if self.os_only:
                for reg in regions:
                    if "_OS" in reg:
                        out += [reg]
            elif self.ss_only:
                for reg in regions:
                    if "_SS" in reg:
                        if reg.startswith("-"):
                            out += [reg[1:]]
                        else:
                            out += ["-" + reg]
            else:
                out = regions
        else:
            for reg in regions:
                anti_sign = "-"
                sign = ""
                reg_nosign = reg
                if reg.startswith("-"):
                    anti_sign = ""
                    sign = "-"
                    reg_nosign = reg[1:]
                if "_OS" in reg_nosign:
                    reg_nosign = reg_nosign.replace("_OS", "_SS")
                else:
                    reg_nosign = reg_nosign.replace("_SS", "_OS")
                if not self.ss_only:
                    out += [reg]
                    out += [f"{anti_sign}{reg_nosign}"]
                else:
                    out += [f"{sign}{reg_nosign}"]
        if not self.loose_sr:
            return out
        else:
            antiSR = []
            for reg in out:
                antiSR += [reg.replace("_SR", "_Anti_SR")]
            if not self.loose_sr_only:
                return antiSR + out
            else:
                return antiSR

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

    def __init__(self, os_ss_sub: bool = False, loose_sr: bool = False, ss_only: bool = False, name: str = ""):
        super().__init__(os_ss_sub=os_ss_sub, loose_sr=loose_sr, ss_only=ss_only)
        if name != "":
            self.name = name

    def get_regions(self, regions):
        return self.format(regions)


class MockMC(ProxyChannel):
    name = 'MockMC'

    def __init__(self, os_ss_sub: bool = False, loose_sr: bool = False, subtract_mj: bool = False):
        super().__init__(os_ss_sub=os_ss_sub, loose_sr=loose_sr)
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

    def __init__(self, os_ss_sub: bool = False, loose_sr: bool = False, fake_factor: bool = False,
                 tight_only: bool = False, loose_only: bool = False):
        super().__init__(os_ss_sub=os_ss_sub, loose_sr=loose_sr)
        self.fake_factor = fake_factor
        self.tight_only = tight_only
        self.loose_only = loose_only
        if self.fake_factor:
            self.name = 'FakeFactor'
        if self.tight_only:
            self.name += '_TightOnly'
        if self.loose_only:
            self.name += '_LooseOnly'

    def get_regions(self, regions):
        out = []
        for reg in regions:
            sign = ""
            anti_sign = "-"
            if reg.startswith("-"):
                reg = reg[1:]
                sign = "-"
                anti_sign = ""
            if not self.tight_only:
                out += [f"{sign}AntiTight_{reg}"]
            if not self.fake_factor and not self.loose_only:
                if not self.tight_only:
                    out += [f"{anti_sign}Tight_{reg}"]
                else:
                    out += [f"{sign}Tight_{reg}"]
        return self.format(out)


class SPGChannel(ProxyChannel):
    name = ""

    def __init__(self, name: str = "", regions_OS: list = [], regions_SS: list = [], always_OS: bool = False):
        super().__init__()
        assert name != ""
        self.name = name
        self.regions_OS = regions_OS
        self.regions_SS = regions_SS
        self.always_OS = always_OS

    def get_regions(self, regions):
        has_OS = False
        has_SS = False
        for reg in regions:
            if "_OS" in reg:
                has_OS = True
            if "_SS" in reg:
                has_SS = True
        if not self.always_OS and (has_OS and not has_SS):
            return self.regions_OS
        elif not self.always_OS and (has_SS and not has_OS):
            return self.regions_SS
        elif not self.always_OS and (has_OS and has_SS):
            return self.regions_OS + [f"-{x}" for x in self.regions_SS]
        elif self.always_OS:
            return self.regions_OS
        else:
            raise Exception("Invalid regions")


class GenericChannel(ProxyChannel):
    name = ""

    def __init__(self, os_ss_sub: bool = False, loose_sr: bool = False, loose_sr_only: bool = False,
                 ss_only: bool = False, os_only: bool = False, force_1tag: bool = False, add_1tag: bool = False, region: str = "", name: str = "",
                 regions_override: list = [], regions_OS: list = [], regions_SS: list = []):
        super().__init__(os_ss_sub=os_ss_sub, loose_sr=loose_sr, loose_sr_only=loose_sr_only,
                         ss_only=ss_only, os_only=os_only, force_1tag=force_1tag, add_1tag=add_1tag)
        self.region = region
        self.regions_override = regions_override
        if name == "":
            self.name = region
        else:
            self.name = name
        self.regions_OS = regions_OS
        self.regions_SS = regions_SS

    def get_regions(self, regions):
        if self.regions_override:
            regions = self.regions_override
        elif self.regions_OS or self.regions_SS:
            has_OS = False
            has_SS = False
            for reg in regions:
                if "_OS" in reg:
                    has_OS = True
                if "_SS" in reg:
                    has_SS = True
            if has_OS and not has_SS:
                regions = self.regions_OS
            elif has_SS and not has_OS:
                regions = self.regions_SS
            elif has_OS and has_SS:
                regions = self.regions_OS + [f"-{x}" for x in self.regions_SS]
        if type(self.region) == str:
            region = f"_{self.region}" if self.region else ""
            return self.format([f"{reg}{region}" for reg in regions])
        else:
            return self.format([f"{reg1}_{reg2}" for reg1 in regions for reg2 in self.region])


class MisMatched(ProxyChannel):
    name = 'MisMatched'
    truth_name = 'MisMatched'

    def __init__(self, os_ss_sub: bool = False, loose_sr: bool = False, pdgId: str = ""):
        super().__init__(os_ss_sub=os_ss_sub, loose_sr=loose_sr)
        self.pdgId = pdgId
        if self.pdgId != "":
            self.name = f"{self.pdgId}{self.name}"
            self.truth_name = self.name
        if loose_sr:
            self.name += "_Loose"

    def get_regions(self, regions):
        return self.format([reg + f"_{self.truth_name}" for reg in regions])


class NoMatchBackground(ProxyChannel):
    name = 'NoMatchBackground'

    def __init__(self, os_ss_sub: bool = False, loose_sr: bool = False, regions_override: bool = False, add_no_truth_match: bool = False):
        super().__init__(os_ss_sub=os_ss_sub, loose_sr=loose_sr)
        self.regions_override = regions_override
        self.add_no_truth_match = add_no_truth_match
        if loose_sr:
            self.name += "_Loose"

    def get_regions(self, regions):
        if self.regions_override:
            regions = self.regions_override
        out = [reg + "_Other" for reg in regions] + [reg + "_HardMisMatched" for reg in regions]
        if self.add_no_truth_match:
            out += regions
        return self.format(out)


class MatchedCharm(ProxyChannel):
    name = 'MatchedCharm'

    def __init__(self, os_ss_sub: bool = False, loose_sr: bool = False, ss_only: bool = False, decayMode: str = "Dplus", name: str = ""):
        super().__init__(os_ss_sub=os_ss_sub, loose_sr=loose_sr, ss_only=ss_only)
        self.decayMode = decayMode
        if name:
            self.name = name
        if loose_sr:
            self.name += "_Loose"

    def get_regions(self, regions):
        if self.decayMode == "Dplus":
            return self.format([reg + "_411MisMatched" for reg in regions] +
                               [reg + "_431MisMatched" for reg in regions] +
                               [reg + "_413MisMatched" for reg in regions] +
                               [reg + "_421MisMatched" for reg in regions] +
                               [reg + "_BaryonMisMatched" for reg in regions])
        elif self.decayMode in ["DstarKPiPi0", "Dstar"]:
            return self.format([reg + "_431MisMatched" for reg in regions] +
                               [reg + "_411MisMatched" for reg in regions] +
                               [reg + "_421MisMatched" for reg in regions] +
                               [reg + "_BaryonMisMatched" for reg in regions])


class MatchedCharmGeom(ProxyChannel):
    name = 'MatchedCharmGeom'

    def __init__(self, os_ss_sub: bool = False, loose_sr: bool = False, ss_only: bool = False, decayMode: str = "Dplus", name: str = ""):
        super().__init__(os_ss_sub=os_ss_sub, loose_sr=loose_sr, ss_only=ss_only)
        self.decayMode = decayMode
        if name:
            self.name = name
        if loose_sr:
            self.name += "_Loose"

    def get_regions(self, regions):
        if self.decayMode == "Dplus":
            return self.format([reg + "_411MisMatchedGeom" for reg in regions] +
                               [reg + "_431MisMatchedGeom" for reg in regions] +
                               [reg + "_413MisMatchedGeom" for reg in regions] +
                               [reg + "_421MisMatchedGeom" for reg in regions] +
                               [reg + "_4132MisMatchedGeom" for reg in regions] +
                               [reg + "_4122MisMatchedGeom" for reg in regions] +
                               [reg + "_4232MisMatchedGeom" for reg in regions] +
                               [reg + "_BaryonMisMatchedGeom" for reg in regions])


class Rest(ProxyChannel):
    name = 'Rest'

    def __init__(self, os_ss_sub: bool = False, loose_sr: bool = False, allowed_regions: list = [],
                 name: str = "", exclude_mismatched: bool = False, include_other: bool = False):
        super().__init__(os_ss_sub=os_ss_sub, loose_sr=loose_sr)
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


class FakeTrack0(ProxyChannel):
    name = "_fakeTrack0"

    def __init__(self, os_ss_sub: bool = False, loose_sr: bool = False, ptbin: int = -1):
        super().__init__(os_ss_sub=os_ss_sub, loose_sr=loose_sr)
        self.ptbin = ptbin
        if self.ptbin >= 0:
            self.name += f"_truth_pt_bin{self.ptbin}"

    def get_regions(self, regions):
        return self.format([reg + f"_{self.name}" for reg in regions])


class FakeTrack1(ProxyChannel):
    name = "_fakeTrack1"

    def __init__(self, os_ss_sub: bool = False, loose_sr: bool = False, ptbin: int = -1):
        super().__init__(os_ss_sub=os_ss_sub, loose_sr=loose_sr)
        self.ptbin = ptbin
        if self.ptbin >= 0:
            self.name += f"_truth_pt_bin{self.ptbin}"

    def get_regions(self, regions):
        return self.format([reg + f"_{self.name}" for reg in regions])


class FakeTrack2(ProxyChannel):
    name = "_fakeTrack2"

    def __init__(self, os_ss_sub: bool = False, loose_sr: bool = False, ptbin: int = -1):
        super().__init__(os_ss_sub=os_ss_sub, loose_sr=loose_sr)
        self.ptbin = ptbin
        if self.ptbin >= 0:
            self.name += f"_truth_pt_bin{self.ptbin}"

    def get_regions(self, regions):
        return self.format([reg + f"_{self.name}" for reg in regions])


class FakeTrack3(ProxyChannel):
    name = "_fakeTrack3"

    def __init__(self, os_ss_sub: bool = False, loose_sr: bool = False, ptbin: int = -1):
        super().__init__(os_ss_sub=os_ss_sub, loose_sr=loose_sr)
        self.ptbin = ptbin
        if self.ptbin >= 0:
            self.name += f"_truth_pt_bin{self.ptbin}"

    def get_regions(self, regions):
        return self.format([reg + f"_{self.name}" for reg in regions])
