#!/usr/bin/env python
import logging
import os
import re
import ROOT
import sys

# ATLAS Style
dirname = os.path.join(os.path.dirname(__file__), "../../atlasrootstyle")
ROOT.gROOT.SetBatch(True)
ROOT.gROOT.LoadMacro(os.path.join(dirname, "AtlasStyle.C"))
ROOT.SetAtlasStyle()

# logging
root = logging.getLogger()
root.setLevel(logging.DEBUG)
handler = logging.StreamHandler(sys.stdout)
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
root.addHandler(handler)

# proxy samples
pdfs = {
    260000: 101,
    261000: 101,
    303600: 101,
    25200: 51,
    65060: 31,
    14400: 59,
    14600: 59,
}

# label dict
label_dict = {
    "meson_eta": "#eta(D)",
    "meson_pt": "p_{T}(D)",
    "lep_eta": "#eta(lep)",
}


def process_file(f, name):
    variations = {}
    at_variation = False
    save_histo = False
    for line in f:
        if "BEGIN HISTO1D" in line:
            logging.info(f"{line}")
            variation = re.findall("\# BEGIN HISTO1D /[\w]+\[(.*)\]", line)  # noqa: W605
            if not len(variation):
                continue
            at_variation = True
            variation = variation[0]
            x_low = []
            x_high = []
            y = []
            y_up = []
            y_dn = []
            print(f"found variation {variation}")
        elif at_variation and "xlow" in line and "xhigh" in line and "val" in line:
            save_histo = True
        elif "END HISTO1D" in line:
            if at_variation:
                variations[variation] = (x_low, x_high, y, y_up, y_dn)
            at_variation = False
            save_histo = False
        elif save_histo:
            vals = [float(x) for x in line.split()]
            x_low += [vals[0]]
            x_high += [vals[1]]
            y += [vals[2]]
            y_up += [vals[3]]
            y_dn += [vals[4]]
    return variations


def main(options):

    for c in options.channels.split(","):

        # make output folder if not exist
        if not os.path.isdir(os.path.join(options.output, c)):
            os.makedirs(os.path.join(options.output, c))

        for v in options.vars.split(","):
            file_name = os.path.join(options.input, f"{c}_{v}.dat")
            f = open(file_name, "r")
            logging.info(f"Succesfully opened file {file_name} as {f}")

            # process the dat file..
            variations = process_file(f, f"{c}_{v}")
            pdf_sets = {}
            for variation in variations:
                if "MUR10" not in variation or "MUF10" not in variation or "DYNSCALE" in variation:
                    continue
                pdf_set = re.findall("PDF([0-9]+)", variation)
                if not len(pdf_set):
                    continue
                pdf_set = int(pdf_set[0])
                for pdf in pdfs:
                    if pdf_set >= pdf and pdf_set <= pdf + pdfs[pdf]:
                        print(f"variation {variation} is part of {pdf}")
                        if pdf not in pdf_sets:
                            pdf_sets[pdf] = {}
                        pdf_sets[pdf][pdf_set] = variations[variation]
            for pdf in pdf_sets:
                logging.info(f"pdf {pdf} has {len(pdf_sets[pdf])} variatins")

            # re-order the data..
            for pdf in pdf_sets:
                logging.info(f"processing pdf {pdf}...")
                f = open(f"{pdf}.txt", "w")
                nominal = pdf_sets[pdf][pdf]
                nbins = len(nominal[2])
                variations = range(pdf, pdf + pdfs[pdf])
                string = 'bin\t' + '\t'.join([str(x) for x in variations]) + '\n'
                f.write(string)
                print(variations)
                for i in range(nbins):
                    values = []
                    for variation in variations:
                        values += [str(pdf_sets[pdf][variation][2][i])]
                    vals = '\t'.join(values)
                    vals = str(i) + '\t' + vals + '\n'
                    f.write(vals)
                f.close()


if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser()

    # ----------------------------------------------------
    # arguments
    # ----------------------------------------------------
    parser.add_option('-i', '--input',
                      action="store", dest="input",
                      help="analysis config file")
    parser.add_option('-o', '--output',
                      action="store", dest="output",
                      help="output folder name",
                      default="rivet_bands")
    parser.add_option('-c', '--channels',
                      action="store", dest="channels",
                      help="run over a subset of channels (comma separated)")
    parser.add_option('-v', '--vars',
                      action="store", dest="vars",
                      help="run over a subset of variables (comma separated)")

    # parse input arguments
    options, args = parser.parse_args()

    # output name
    out_name = options.output

    # make output folder if not exist
    if not os.path.isdir(out_name):
        os.makedirs(out_name)

    # do the plotting
    main(options)
