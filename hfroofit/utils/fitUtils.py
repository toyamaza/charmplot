import ROOT
import json
import os

# C++ implementations
cpp_name = os.path.join(os.path.dirname(__file__), "fitUtils.cxx")
ROOT.gROOT.LoadMacro(cpp_name)


def run_fit(w, datasetName):
    # load model from workspace
    mc = w.genobj("ModelConfig")
    pdf = w.pdf(mc.GetPdf().GetName())

    # load data from workspace
    data = w.data(datasetName)

    # do the fit
    res = pdf.fitTo(data, ROOT.RooFit.SumW2Error(True), ROOT.RooFit.Save(True))
    return res


def results_to_dict(res, pars):
    out = dict()
    for par in pars:
        p = res.floatParsFinal().find(par)
        if not p:
            continue
        out.update({par: [p.getVal(), p.getAsymErrorHi()]})
    return out


def dict_to_json(dict, path, name):
    pfn = os.path.join(os.path.dirname(path), name + ".json")
    with open(pfn, 'w') as outfile:
        json.dump(dict, outfile)
