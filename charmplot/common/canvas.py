from charmplot.common import utils
from charmplot.common import likelihoodFit
from charmplot.control import channel
from charmplot.control import variable
from charmplot.control import sample
from typing import Dict
import logging
import math
import ROOT
import sys

MC_Map = Dict[sample.Sample, ROOT.TH1]

# logging
logger = logging.getLogger(__name__)

BASE_WIDTH = 800.
BASE_HEIGHT = 600.


class CanvasBase(object):

    def __init__(self, c: channel.Channel, v: variable.Variable, x: float, y: float):
        self.name = c.name
        self.channel = c
        self.variable = v
        self.x = x
        self.y = y
        self.canv = ROOT.TCanvas(self.name, "", int(self.x), int(self.y))
        self.set_canvas_margins()

        # ATLAS label and text
        self.text_height = 32. / self.y
        self.text_height_small = 28. / self.y
        self.text_pos_y = 1 - 2.5 * self.text_height
        logger.debug(f"created canvas with size {self.x} {self.y}")

    def set_x_range(self, h):
        if self.variable.x_range:
            h.GetXaxis().SetRangeUser(self.variable.x_range[0], self.variable.x_range[1])

    def make_proxy_histogram(self, h, name="proxy"):
        proxy = h.Clone(f"{self.name}_{name}")
        for i in range(0, proxy.GetNbinsX() + 2):
            proxy.SetBinContent(i, 0)
            proxy.SetBinError(i, 0)
        proxy.Draw("hist")
        return proxy

    def set_canvas_margins(self):
        self.canv.SetTopMargin(40 / self.canv.GetWindowHeight())
        self.canv.SetBottomMargin(80 / self.canv.GetWindowHeight())
        self.canv.SetLeftMargin(120 / self.canv.GetWindowWidth())
        self.canv.SetRightMargin(30 / self.canv.GetWindowWidth())

    def set_axis_title(self, proxy, title="", events="Events"):
        self.bin_width = proxy.GetBinWidth(1)
        precision = 1
        bin_width = self.bin_width
        while (bin_width % 1 > 1e-4 and precision < 4):
            precision += 1
            bin_width *= 10

        if self.variable.label:
            if self.variable.unit:
                proxy.GetXaxis().SetTitle(f"{self.variable.label} [{self.variable.unit}]")
                proxy.GetYaxis().SetTitle(f"{events} / ({self.bin_width:.{precision}f} {self.variable.unit})")
            else:
                proxy.GetXaxis().SetTitle(f"{self.variable.label}")
                proxy.GetYaxis().SetTitle(f"{events} / {self.bin_width:.{precision}f}")
        else:
            proxy.GetXaxis().SetTitle(f"{self.variable.name}")
            proxy.GetYaxis().SetTitle(f"{events} / {self.bin_width:.{precision}f}")

        if title:
            proxy.GetYaxis().SetTitle(title)

    def text(self, text):
        line = ROOT.TLatex()
        line.SetTextFont(43)
        line.SetTextSize(32)
        line.DrawLatex(0.18, self.text_pos_y, text)
        self.text_pos_y -= self.text_height

    def set_axis_text_size(self, proxy, y=1.0, no_x_axis=False):
        # y-axis
        proxy.GetYaxis().SetTitleFont(43)
        proxy.GetYaxis().SetTitleSize(32)
        proxy.GetYaxis().SetTitleOffset(1.4 * (self.y / BASE_HEIGHT))
        proxy.GetYaxis().SetLabelFont(43)
        proxy.GetYaxis().SetLabelSize(32)
        proxy.GetYaxis().SetLabelOffset(0.005)

        # x-axis
        proxy.GetXaxis().SetTitleFont(43)
        proxy.GetXaxis().SetTitleSize(32)
        proxy.GetXaxis().SetTitleOffset(1.0 / y)
        proxy.GetXaxis().SetLabelFont(43)
        proxy.GetXaxis().SetLabelSize(32)
        proxy.GetXaxis().SetLabelOffset(0.005)
        if no_x_axis:
            proxy.GetXaxis().SetTitleSize(0)
            proxy.GetXaxis().SetLabelSize(0)

    def print(self, name):
        self.canv.Print(name)


class Canvas2(CanvasBase):

    fit = None
    pad1 = None
    pad2 = None
    legend = None
    offset = 0

    def __init__(self, c: channel.Channel, v: variable.Variable, x: float, y: float, y_split: float = 0.35, fit: likelihoodFit.LikelihoodFit = None):
        super(Canvas2, self).__init__(c, v, x, y)

        # upper/lower canvas split
        self.y_split = y_split
        if self.y_split:
            self.offset = 0.05

        # fit results
        self.fit = fit

        # upper pad
        pad1 = ROOT.TPad(self.name + "_upper_pad", "", 0., self.y_split, 1., 1.)
        pad1.SetFillStyle(4000)
        pad1.SetTopMargin(self.canv.GetTopMargin() / (1 - self.y_split))
        pad1.SetLeftMargin(self.canv.GetLeftMargin())
        pad1.SetRightMargin(self.canv.GetRightMargin())
        if self.y_split:
            pad1.SetBottomMargin(self.offset)
        else:
            pad1.SetBottomMargin(self.canv.GetBottomMargin())
        pad1.Draw()
        self.pad1 = pad1

        # set ATLAS label
        self.set_atlas_label()

        # print fit results
        self.fit_results()

        # lower pad
        if self.y_split:
            pad2 = ROOT.TPad(self.name + "_upper_pad", "", 0., 0, 1., 1 - (1 - self.y_split) * (1 - self.offset))
            pad2.SetGridy()
            pad2.SetFillStyle(4000)
            pad2.SetTopMargin(0)
            pad2.SetBottomMargin(self.canv.GetBottomMargin() / (self.y_split + self.offset))
            pad2.SetLeftMargin(self.canv.GetLeftMargin())
            pad2.SetRightMargin(self.canv.GetRightMargin())
            pad2.Draw()
            self.pad2 = pad2

    def set_atlas_label(self):
        # ATLAS label
        self.atlas_label("internal")
        self.text("#sqrt{s} = 13 TeV, %.1f fb^{-1}" % (utils.get_lumi(self.channel.lumi) / 1000.))
        for label in self.channel.label:
            self.text(f"{label}")

    def fit_results(self):
        if self.fit:
            for s, r in self.fit.result.items():
                if 'Multijet' in s or s in self.fit.fixed:
                    continue
                self.text(f"#mu({s}) = {r[0]:.3f}")

    def print_all(self, output, channel, var, multipage_pdf=False, first_plot=False, last_plot=False, as_png=False):
        self.pad1.cd()
        ROOT.gPad.RedrawAxis()
        self.pad2.cd()
        ROOT.gPad.RedrawAxis()
        self.print(f"{output}/{channel}/{channel}_{var}.pdf")
        if as_png:
            self.print(f"{output}/{channel}/{channel}_{var}.png")
        if multipage_pdf:
            if first_plot:
                self.print(f"{output}/{channel}.pdf(")
            else:
                self.print(f"{output}/{channel}.pdf")
        self.set_logy()
        self.print(f"{output}/{channel}/{channel}_{var}_LOG.pdf")
        if as_png:
            self.print(f"{output}/{channel}/{channel}_{var}_LOG.png")
        if multipage_pdf:
            if last_plot:
                self.print(f"{output}/{channel}.pdf)")
            else:
                self.print(f"{output}/{channel}.pdf")

    def configure_histograms(self, mc_map: MC_Map, data: ROOT.TH1, v: variable.Variable):
        if data:
            data.SetMarkerSize(0.8)
        for s, h in mc_map.items():
            logger.debug(f"configuring {s} {h}")
            if s.fillColor:
                h.SetFillColor(s.fillColor)
            else:
                h.SetFillStyle(0)
            h.SetLineWidth(0)

    def set_maximum(self, histograms: list, variable: variable.Variable, mc_min: ROOT.TH1 = None):

        if not histograms:
            logger.critical("No input histograms given.")
            sys.exit(1)

        # x range
        x_min = histograms[0].GetBinCenter(1)
        x_max = histograms[0].GetBinCenter(histograms[0].GetNbinsX())
        if variable.x_range:
            x_min = variable.x_range[0]
            x_max = variable.x_range[1]

        # Get maximum on both sides of the plot
        max_left_ = 0
        max_right_ = 0
        for h in histograms:
            max_left = utils.get_maximum(h, x_min, x_min + 0.75 * (x_max - x_min))
            max_right = utils.get_maximum(h, x_min + 0.75 * (x_max - x_min), x_max - h.GetBinWidth(h.GetNbinsX()))
            if max_left > max_left_:
                max_left_ = max_left
            if max_right > max_right_:
                max_right_ = max_right

        # Specific minimum value based on some mc histogram
        self.max_val = max(max_left_, max_right_)
        if mc_min:
            self.min_positive_val = self.max_val
            for i in range(1, mc_min.GetNbinsX() + 1):
                y = mc_min.GetBinContent(i)
                if y > 0 and y < self.min_positive_val:
                    self.min_positive_val = y

        # max y on left side
        max_left_y = ((self.text_pos_y - self.y_split) / (1 - self.y_split - self.pad1.GetTopMargin()))

        # max y on right side
        max_right_y = None
        if self.legend:
            max_right_y = self.leg_y1

        # Determine whether maximum is on the left or the right side of the plot
        if not max_right_y or max_right_ <= 0 or (max_left_ > 0 and max_left_ / max_right_ > 1.5):
            self.maximum_scale_factor = 1.05 / max_left_y
        else:
            self.maximum_scale_factor = 1.05 / min(max_right_y, max_left_y)
        self.proxy_up.SetMaximum(self.maximum_scale_factor * self.max_val)

    def set_logy(self):
        if not hasattr(self, "max_val"):
            logger.critical("please call 'set_maximum' before setting logy")
            sys.exit(1)

        if hasattr(self, "min_positive_val"):
            self.proxy_up.SetMinimum(0.5 * self.min_positive_val)
        self.pad1.cd()
        self.pad1.SetLogy()
        self.proxy_up.SetMaximum(math.pow(10, math.log10(self.max_val) * self.maximum_scale_factor +
                                          (1 - self.maximum_scale_factor) * math.log10(self.proxy_up.GetMinimum())))

    def make_legend(self, data, mc_tot, mc_map, samples, print_yields=False):
        # temp entry for sys unc
        temp_err = ROOT.TGraphErrors()
        temp_err.SetLineColor(ROOT.kBlack)
        temp_err.SetFillColor(ROOT.kGray + 2)
        temp_err.SetFillStyle(3354)
        self.temp_err = temp_err

        # legend
        n_entries = len(samples)
        if data:
            n_entries += 1
        if mc_tot:
            n_entries += 1
        self.leg_y2 = 1 - 1.8 * self.text_height_small / (1 - self.y_split)
        self.leg_y1 = self.leg_y2 - n_entries * self.text_height_small / (1 - self.y_split)
        leg = ROOT.TLegend(0.65, self.leg_y1, 0.9, self.leg_y2)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(28)
        leg.SetTextFont(43)
        if data:
            if print_yields:
                leg.AddEntry(data, "Data #scale[0.50]{#splitline{%.2e}{/ MC = %1.3f}}" % (data.GetSum(), data.GetSum() / mc_tot.GetSum()), "pe")
            else:
                leg.AddEntry(data, "Data", "pe")
        if mc_tot:
            if print_yields:
                leg.AddEntry(temp_err, "SM tot. #scale[0.50]{%.2e}" % mc_tot.GetSum(), "lf")
            else:
                leg.AddEntry(temp_err, "SM tot.", "lf")
        for s in samples:
            name = s.name
            if hasattr(s, "legendLabel"):
                name = s.legendLabel
            if print_yields:
                leg.AddEntry(mc_map[s], "%s #scale[0.50]{%.2e}" % (name, mc_map[s].GetSum()), "f")
            else:
                leg.AddEntry(mc_map[s], "%s" % name, "f")
        self.pad1.cd()
        self.legend = leg
        self.legend.Draw()

    def atlas_label(self, text):
        l1 = ROOT.TLatex()
        l1.SetTextFont(73)
        l1.SetTextSize(32)
        l1.DrawLatex(0.18, self.text_pos_y, "ATLAS")
        if text:
            l2 = ROOT.TLatex()
            l2.SetTextFont(43)
            l2.SetTextSize(32)
            l2.DrawLatex(0.18 + 0.14, self.text_pos_y, text)
        self.text_pos_y -= self.text_height

    def set_ratio_range(self, dn, up, ndivisions=506):
        if self.variable.ratio_range:
            self.proxy_dn.SetMinimum(self.variable.ratio_range[0])
            self.proxy_dn.SetMaximum(self.variable.ratio_range[1])
        else:
            self.proxy_dn.SetMinimum(dn)
            self.proxy_dn.SetMaximum(up)
        self.proxy_dn.GetYaxis().SetNdivisions(ndivisions)

    def draw_qcd_frac(self, h, gr_err):
        self.pad2.cd()
        h.Draw("same hist")
        gr_err.Draw("same e2")
        self.proxy_dn.SetMinimum(0)
        self.proxy_dn.SetMaximum(0.99)
        self.proxy_dn.GetYaxis().SetTitle("QCD fraction")

    def construct(self, h):
        # proxy histogram to control the upper axis
        self.pad1.cd()
        self.proxy_up = self.make_proxy_histogram(h, "up")
        self.set_axis_title(self.proxy_up)
        self.set_axis_text_size(self.proxy_up, no_x_axis=(self.y_split > 0))
        self.set_x_range(self.proxy_up)

        # proxy histogram to control the lower axis
        if self.pad2:
            self.pad2.cd()
            self.proxy_dn = self.make_proxy_histogram(h, "dn")
            self.set_axis_title(self.proxy_dn, "Data / MC" if self.y_split else "")
            self.set_axis_text_size(self.proxy_dn, self.y_split + self.offset)
            self.set_ratio_range(0.75, 1.24)
            self.set_x_range(self.proxy_dn)


class CanvasMCRatio(Canvas2):

    def __init__(self, c: channel.Channel, v: variable.Variable, x: float, y: float, y_split: float = 0.35):
        super(CanvasMCRatio, self).__init__(c, v, x, y, y_split)

    def make_legend(self, mc_map, samples):
        n_entries = len(samples)
        self.leg_y2 = 1 - 1.8 * self.text_height_small / (1 - self.y_split)
        self.leg_y1 = self.leg_y2 - n_entries * self.text_height_small / (1 - self.y_split)
        leg = ROOT.TLegend(0.65, self.leg_y1, 0.9, self.leg_y2)
        leg.SetBorderSize(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(28)
        leg.SetTextFont(43)
        for s in samples:
            name = s.name
            if hasattr(s, "legendLabel"):
                name = s.legendLabel
            leg.AddEntry(mc_map[s], "%s" % name, "l")
        self.pad1.cd()
        self.legend = leg
        self.legend.Draw()

    def configure_histograms(self, mc_map: MC_Map, v: variable.Variable):
        for s, h in mc_map.items():
            logger.debug(f"configuring {s} {h}")
            h.SetFillStyle(0)
            if s.lineColor:
                h.SetLineColor(s.lineColor)
                h.SetLineWidth(2)
            h.Scale(1. / h.GetSum())

    def set_atlas_label(self):
        # ATLAS label
        self.atlas_label("internal")
        self.text("#sqrt{s} = 13 TeV")
        for label in self.channel.label:
            self.text(f"{label}")

    def construct(self, h):
        # proxy histogram to control the upper axis
        self.pad1.cd()
        self.proxy_up = self.make_proxy_histogram(h, "up")
        self.set_axis_title(self.proxy_up, events="Normalized Entries")
        self.set_axis_text_size(self.proxy_up, no_x_axis=(self.y_split > 0))
        self.set_x_range(self.proxy_up)

        # proxy histogram to control the lower axis
        if self.pad2:
            self.pad2.cd()
            self.proxy_dn = self.make_proxy_histogram(h, "dn")
            self.set_axis_title(self.proxy_dn, "Data / MC" if self.y_split else "")
            self.set_axis_text_size(self.proxy_dn, self.y_split + self.offset)
            self.set_ratio_range(0.50, 1.49)
            self.set_x_range(self.proxy_dn)
