#!/usr/bin/env python
from copy import deepcopy
from math import log10, floor
import os
import ROOT
import yaml

ROOT.gROOT.SetBatch(True)

DPLUS_FOLDER = "/global/cscratch1/sd/mmuskinj/TRExFitter/fit_2022_12_14_Dplus_new/"
DSTAR_FOLDER = "/global/cscratch1/sd/mmuskinj/TRExFitter/fit_2022_12_14_Dstar_new/"
# DPLUS_FOLDER = "/global/cfs/cdirs/atlas/wcharm/TRExFitter/Output/Dplus_2022_08_05_v2"
# DSTAR_FOLDER = "/global/cfs/cdirs/atlas/wcharm/TRExFitter/Output/Dstar_2022_08_11_v2"

POIs_abs = [f"mu_Wplus_{i}" for i in range(1, 6)] + [f"mu_Wminus_{i}" for i in range(1, 6)]

PRUNING_THRESHOLD = 0.98

QUALIFIERS = [{'name': 'SQRT(s)', 'units': 'GeV', 'value': 13000}, {'name': 'LUMINOSITY', 'units': 'fb$^{-1}$', 'value': 140.1}]

CHANNELS = ["dplus", "dstar"]

NP_DESCRIPTION = {
    "EG_RESOLUTION_ALL": "Electron Resolution",
    "EG_SCALE_ALL": "Electron Energy Resolution",
    "EL_CHARGEID_STAT": "Electron Charge ID (Stat.)",
    "EL_CHARGEID_SYStotal": "Electron Charge ID (Syst.)",
    "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR": "Electron ID Efficiency",
    "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR": "Electron Isolation Efficiency",
    "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR": "Electron Reco Efficiency",
    "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR": "Electron Trigger Efficiency",
    "FID_EFF_AS": "Fiducial Efficiency AlphaS Variation",
    "FID_EFF_PDF": "Fiducial Efficiency PDF Variation",
    "FID_EFF_QCD": "Fiducial Efficiency Scale Variation",
    "FID_EFF_WPLUSD_MG": "Fiducial Efficiency MC Variation",
    "JET_EffectiveNP_8restTerm": "Jet Energy Scale NP 8 (rest term)",
    "JET_EtaIntercalibration_Modelling": "Jet Intercalibration Modelling",
    "JET_EtaIntercalibration_NonClosure_highE": "Jet Intercalibration Closure (high E)",
    "JET_EtaIntercalibration_NonClosure_negEta": "Jet Intercalibration Closure (neg eta)",
    "JET_EtaIntercalibration_NonClosure_posEta": "Jet Intercalibration Closure (pos eta)",
    "JET_EtaIntercalibration_TotalStat": "Jet Intercalibration Stat.",
    "JET_Flavor_Composition": "Jet Flavor Composition",
    "JET_Flavor_Response": "Jet Flavor Response",
    "JET_JER_DataVsMC_MC16": "Jet Energy Resolution Data vs MC",
    "JET_JER_EffectiveNP_7restTerm": "Jet Energy Resolution NP 7 (rest term)",
    "JET_JvtEfficiency": "Jet JVT Efficiency",
    "JET_Pileup_OffsetMu": "Jet Energy Scale PU Offset (Mu)",
    "JET_Pileup_OffsetNPV": "Jet Energy Scale PU Offset (NPV)",
    "JET_Pileup_PtTerm": "Jet Energy Scale PU pT Term",
    "JET_Pileup_RhoTopology": "Jet Pileup Rho Topology",
    "JET_PunchThrough_MC16": "Jet Flavor Response",
    "Lumi": "Luminosity",
    "MET_SoftTrk_ResoPara": "MET Soft Term Parallel Resolution",
    "MET_SoftTrk_ResoPerp": "MET Soft Term Perpendicular Resolution",
    "MET_SoftTrk": "MET Soft Term Scale",
    "mu_Top": "Normalization Factor for Top",
    "MUON_CB": "Muon Energy Resolution",
    "MUON_EFF_ISO_STAT": "Muon Isolation Efficiency (Stat.)",
    "MUON_EFF_ISO_SYS": "Muon Isolation Efficiency (Syst.)",
    "MUON_EFF_RECO_STAT_LOWPT": "Muon Reco Eff. Low pT (Stat.)",
    "MUON_EFF_RECO_STAT": "Muon Reco Efficiency (Stat.)",
    "MUON_EFF_RECO_SYS_LOWPT": "Muon Reco Eff. Low pT (Syst.)",
    "MUON_EFF_RECO_SYS": "Muon Reco Efficiency (Syst.)",
    "MUON_EFF_TrigStatUncertainty": "Muon Trigger Efficiency (Stat].)",
    "MUON_EFF_TrigSystUncertainty": "Muon Trigger Efficiency (Syst.)",
    "MUON_EFF_TTVA_STAT": "Muon TTVA Efficiency (Stat.)",
    "MUON_EFF_TTVA_SYS": "Muon TTVA Efficiency (Syst.)",
    "MUON_SAGITTA_DATASTAT": "Muon Sagitta Bias (Stat.)",
    "MUON_SAGITTA_RESBIAS": "Muon Sagitta Bias (Syst.)",
    "MUON_SCALE": "Muon Energy Scale",
    "PROD_FRAC_EIG_1": "Charm Production Fraction NP 1",
    "PROD_FRAC_EIG_2": "Charm Production Fraction NP 2",
    "PROD_FRAC_EIG_3": "Charm Production Fraction NP 3",
    "PRW_DATASF": "Pileup Reweighting",
    "Top_GEN_isr": "Top Bkg ISR Uncertainty",
    "Top_GEN_muF": "Top Bkg muF Uncertainty",
    "Top_GEN_muR": "Top Bkg muR Uncertainty",
    "Top_GEN_PDF": "Top Bkg PDF Uncertainty",
    "Top_GEN_Var3c": "Top Bkg Var3c Variation",
    "Top_HS": "Top Bkg Hard Scatter Variation",
    "Top_Shower": "Top Bkg Hard Shower Variation",
    "TRK_EFF_D0_DEAD": "Track d0 Resolution (Dead Pixels)",
    "TRK_EFF_D0_MEAS": "Track d0 Resolution",
    "TRK_EFF_IBL": "IBL Track Efficiency",
    "TRK_EFF_Overal": "Overall Tracking Efficiency",
    "TRK_EFF_PP0": "PP0 Track Efficiency",
    "TRK_EFF_QGSP": "Physics List Track Efficiency",
    "TRK_EFF_Z0_DEAD": "Track z0 Resolution (Dead Pixels)",
    "TRK_EFF_Z0_MEAS": "Track z0 Resolution",
    "TRK_EFF_D0_DEAD_Top": "Track d0 Resolution (Dead Pixels; Top)",
    "TRK_EFF_D0_MEAS_Top": "Track d0 Resolution (Top)",
    "TRK_EFF_IBL_Top": "IBL Track Efficiency (Top)",
    "TRK_EFF_Overal_Top": "Overall Tracking Efficiency (Top)",
    "TRK_EFF_PP0_Top": "PP0 Track Efficiency (Top)",
    "TRK_EFF_QGSP_Top": "Physics List Track Efficiency (Top)",
    "TRK_EFF_Z0_DEAD_Top": "Track z0 Resolution (Dead Pixels; Top)",
    "TRK_EFF_Z0_MEAS_Top": "Track z0 Resolution (Top)",
}
NP_DESCRIPTION.update({f"FT_EFF_Eigen_B_{i}": f"FTAG Efficinecy (b-jet NP {i+1})" for i in range(0, 9)})
NP_DESCRIPTION.update({f"FT_EFF_Eigen_C_{i}": f"FTAG Efficinecy (c-jet NP {i+1})" for i in range(0, 9)})
NP_DESCRIPTION.update({f"FT_EFF_Eigen_Light_{i}": f"FTAG Efficinecy (u-jet NP {i+1})" for i in range(0, 9)})
NP_DESCRIPTION.update({f"JET_EffectiveNP_{i}": f"Jet Energy Scale NP {i}" for i in range(1, 10)})
NP_DESCRIPTION.update({f"JET_JER_EffectiveNP_{i}": f"Jet Energy Resolution NP {i}" for i in range(1, 7)})
NP_DESCRIPTION.update({
    "DPLUS_BR_BKG": "D+ to non-Kpipi Branching Ratio Unc.",
    "DPLUS_BR_SIG": "D+ to Kpipi Branching Ratio Unc.",
    "DPLUS_Matched_Sherpa_MG_2P": "W+D Shape Unc. (D+)",
    "DPLUS_MM_EL_FR__Wjets_MadGraph": "MM El. Fake Rate Modelling (D+)",
    "DPLUS_MM_EL_FR_MET_Var": "MM El. Fake Rate MET Variation (D+)",
    "DPLUS_MM_EL_FR_STAT": "MM El. Fake Rate Stat. (D+)",
    "DPLUS_MM_EL_RR_MadGraph_Var": "MM El. Real Rate Modelling (D+)",
    "DPLUS_MM_EL_RR_STAT": "MM El. Real Rate Stat. (D+)",
    "DPLUS_MM_MU_FR__Wjets_MadGraph": "MM Mu. Fake Rate Modelling (D+)",
    "DPLUS_MM_MU_FR_MET_Var": "MM Mu. Fake Rate MET Variation (D+)",
    "DPLUS_MM_MU_FR_STAT": "MM Mu. Fake Rate Stat. (D+)",
    "DPLUS_MM_MU_RR_MadGraph_Var": "MM Mu. Real Rate Modelling (D+)",
    "DPLUS_MM_MU_RR_STAT": "MM Mu. Real Rate Stat. (D+)",
    "DPLUS_PEAK_POSITION_MINUS": "D- Peak Position",
    "DPLUS_PEAK_POSITION_PLUS": "D+ Peak Position",
    "DPLUS_PEAK_RESOLUTION_MINUS": "D- Peak Resolution",
    "DPLUS_PEAK_RESOLUTION_PLUS": "D+ Peak Resolution",
    "DPLUS_SHERPA2211_AS_WCHARM": "W+c Bkg AlphaS Uncertainty (D+)",
    "DPLUS_SHERPA2211_AS_WJETS": "W+jets Bkg AlphaS Uncertainty (D+)",
    "DPLUS_SHERPA2211_PDF_WCHARM": "W+c Bkg PDF Uncertainty (D+)",
    "DPLUS_SHERPA2211_PDF_WJETS": "W+jets Bkg PDF Uncertainty (D+)",
    "DPLUS_SHERPA2211_QCD_WCHARM": "W+c Bkg Scale Uncertainty (D+)",
    "DPLUS_SHERPA2211_QCD_WJETS": "W+jets Bkg Scale Uncertainty (D+)",
    "DSTAR_BR_SIG": "D* to (Kpi)pi Branching Ratio Unc.",
    "DSTAR_MM_EL_FR__Wjets_MadGraph": "MM El. Fake Rate Modelling (D*)",
    "DSTAR_MM_EL_FR_MET_Var": "MM El. Fake Rate MET Variation (D*)",
    "DSTAR_MM_EL_FR_STAT": "MM El. Fake Rate Stat. (D*)",
    "DSTAR_MM_EL_RR_MadGraph_Var": "MM El. Real Rate Modelling (D*)",
    "DSTAR_MM_EL_RR_STAT": "MM El. Real Rate Stat. (D*)",
    "DSTAR_MM_MU_FR__Wjets_MadGraph": "MM Mu. Fake Rate Modelling (D*)",
    "DSTAR_MM_MU_FR_MET_Var": "MM Mu. Fake Rate MET Variation (D*)",
    "DSTAR_MM_MU_FR_STAT": "MM Mu. Fake Rate Stat. (D*)",
    "DSTAR_MM_MU_RR_MadGraph_Var": "MM Mu. Real Rate Modelling (D*)",
    "DSTAR_MM_MU_RR_STAT": "MM Mu. Real Rate Stat. (D*)",
    "DSTAR_PEAK_POSITION_ADHOC_MINUS": "D*- Peak Position (0.1 MeV)",
    "DSTAR_PEAK_POSITION_ADHOC_PLUS": "D*+ Peak Position (0.1 MeV)",
    "DSTAR_PEAK_POSITION_MINUS": "D*- Peak Position",
    "DSTAR_PEAK_POSITION_PLUS": "D*+ Peak Position",
    "DSTAR_PEAK_RESOLUTION_MINUS": "D*- Peak Resolution",
    "DSTAR_PEAK_RESOLUTION_PLUS": "D*+ Peak Resolution",
    "DSTAR_SHERPA2211_AS_WJETS": "W+jets Bkg AlphaS Uncertainty (D*)",
    "DSTAR_SHERPA2211_AS_WPLUSC": "W+c Bkg AlphaS Uncertainty (D*)",
    "DSTAR_SHERPA2211_PDF_WJETS": "W+jets Bkg PDF Uncertainty (D*)",
    "DSTAR_SHERPA2211_PDF_WPLUSC": "W+c Bkg PDF Uncertainty (D*)",
    "DSTAR_SHERPA2211_QCD_WJETS": "W+jets Bkg Scale Uncertainty (D*)",
    "DSTAR_SHERPA2211_QCD_WPLUSC": "W+c Bkg Scale Uncertainty (D*)",
    "PRW_DATASF_WJETS": "Pileup Reweighting (W+jets)",
    "WJETS_FIT": "W+jets Parametric Fit Stat. (D*)",
})
NP_DESCRIPTION.update({"DPLUS_CharmNorm_1tag": "W+c(match) Norm. Unc. (D+ 1-tag)"})
NP_DESCRIPTION.update({"DPLUS_MisMatchNorm_1tag": "W+c(mis-match) Norm. Unc. (D+ 1-tag)"})
NP_DESCRIPTION.update({"DPLUS_OtherNorm_0tag": "Other Norm. Unc. (D+ 0-tag)"})
NP_DESCRIPTION.update({"DPLUS_OtherNorm_1tag": "Other Norm. Unc. (D+ 1-tag)"})
NP_DESCRIPTION.update({"DPLUS_WjetsNorm_1tag": "W+jets Norm. Unc. (D+ 1-tag)"})
NP_DESCRIPTION.update({f"DPLUS_CharmNorm_0tag_bin{i}": f"W+c(match) Norm. Unc. (D+ bin {i})" for i in range(1, 6)})
NP_DESCRIPTION.update({f"DPLUS_MisMatch_Sh_MG_2P_bin{i}": f"W+c(mis-match) Shape Unc. (D+ bin {i})" for i in range(1, 6)})
NP_DESCRIPTION.update({f"DPLUS_MisMatchNorm_0tag_bin{i}": f"W+c(mis-match) Norm Unc. (D+ bin {i})" for i in range(1, 6)})
NP_DESCRIPTION.update({f"DPLUS_Wjets_Sh_MG_2P_bin{i}": f"W+jets Shape Unc. (D+ bin {i})" for i in range(1, 6)})
NP_DESCRIPTION.update({f"DPLUS_WjetsNorm_0tag_bin{i}": f"W+jets Norm. Unc. (D+ bin {i})" for i in range(1, 6)})
NP_DESCRIPTION.update({"DSTAR_CharmNorm_1tag": "W+c(match) Norm. Unc. (D* 1-tag)"})
NP_DESCRIPTION.update({"DSTAR_MisMatchNorm_1tag": "W+c(mis-match) Norm. Unc. (D* 1-tag)"})
NP_DESCRIPTION.update({"DSTAR_OtherNorm_0tag": "Other Norm. Unc. (D* 0-tag)"})
NP_DESCRIPTION.update({"DSTAR_OtherNorm_1tag": "Other Norm. Unc. (D* 1-tag)"})
NP_DESCRIPTION.update({"DSTAR_WjetsNorm_1tag": "W+jets Norm. Unc. (D* 1-tag)"})
NP_DESCRIPTION.update({f"DSTAR_CharmNorm_0tag_bin{i}": f"W+c(match) Norm. Unc. (D* bin {i})" for i in range(1, 6)})
NP_DESCRIPTION.update({f"DSTAR_MisMatched_Sherpa_MG_2P_bin{i}": f"W+c(mis-match) Shape. Unc. (D* bin {i})" for i in range(1, 6)})
NP_DESCRIPTION.update({f"DSTAR_MisMatchNorm_0tag_bin{i}": f"W+c(mis-match) Norm Unc. (D* bin {i})" for i in range(1, 6)})
NP_DESCRIPTION.update({f"DSTAR_Wjets_MG_Sherpa_2P_bin{i}": f"W+jets Shape. Unc. (D* bin {i})" for i in range(1, 6)})
NP_DESCRIPTION.update({f"DSTAR_WjetsNorm_0tag_bin{i}": f"W+jets Norm Unc. (D* bin {i})" for i in range(1, 6)})
NP_DESCRIPTION.update({f"DPLUS_Stat_{charge}_DibosonZjets_0tag_{OSSS}_{var}_bin{diff_bin}": f"Other Norm. Stat. D+ (W{('-' if charge == 'minus' else '+')} {var}{diff_bin})" for charge in [
    "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_DESCRIPTION.update({f"DPLUS_Stat_{charge}_WCharm_0tag_{OSSS}_{var}_bin{diff_bin}": f"W+c(match) Norm. Stat. D+ (W{('-' if charge == 'minus' else '+')} {var}{diff_bin})" for charge in [
    "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_DESCRIPTION.update({f"DPLUS_Stat_{charge}_Wjets_0tag_{OSSS}_{var}_bin{diff_bin}": f"W+jets Norm. Stat. D+ (W{('-' if charge == 'minus' else '+')} {var}{diff_bin})" for charge in [
    "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_DESCRIPTION.update({f"DPLUS_Stat_{charge}_WMisMatched_0tag_{OSSS}_{var}_bin{diff_bin}": f"W+c(mis-match) Norm. Stat. D+ (W{('-' if charge == 'minus' else '+')} {var}{diff_bin})" for charge in [
    "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_DESCRIPTION.update({f"DPLUS_Stat_{charge}_Top_0tag_{OSSS}_{var}_bin{diff_bin}": f"Top Norm. Stat. D+ (W{('-' if charge == 'minus' else '+')} {var}{diff_bin})" for charge in [
    "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_DESCRIPTION.update({f"DSTAR_Stat_{charge}_DibosonZjets_0tag_{OSSS}_{var}_bin{diff_bin}": f"Other Norm. Stat. D+ (W{('-' if charge == 'minus' else '+')} {var}{diff_bin})" for charge in [
    "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_DESCRIPTION.update({f"DSTAR_Stat_{charge}_WCharm_0tag_{OSSS}_{var}_bin{diff_bin}": f"W+c(match) Norm. Stat. D+ (W{('-' if charge == 'minus' else '+')} {var}{diff_bin})" for charge in [
    "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_DESCRIPTION.update({f"DSTAR_Stat_{charge}_Wjets_0tag_{OSSS}_{var}_bin{diff_bin}": f"W+jets Norm. Stat. D+ (W{('-' if charge == 'minus' else '+')} {var}{diff_bin})" for charge in [
    "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_DESCRIPTION.update({f"DSTAR_Stat_{charge}_WMisMatched_0tag_{OSSS}_{var}_bin{diff_bin}": f"W+c(mis-match) Norm. Stat. D+ (W{('-' if charge == 'minus' else '+')} {var}{diff_bin})" for charge in [
    "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_DESCRIPTION.update({f"DSTAR_Stat_{charge}_Top_0tag_{OSSS}_{var}_bin{diff_bin}": f"Top Norm. Stat. D+ (W{('-' if charge == 'minus' else '+')} {var}{diff_bin})" for charge in [
    "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_DESCRIPTION.update({f"DPLUS_MC_Stat_Fid_Eff_{charge}_{i}_{j}": f"Response Matrix Stat. D+ (W{('-' if charge == 'minus' else '+')} {i}, {j})" for i in range(1, 6)
                       for j in range(1, 6) for charge in ["minus", "plus"]})
NP_DESCRIPTION.update({f"DSTAR_MC_Stat_Fid_Eff_{charge}_{i}_{j}": f"Response Matrix Stat. D* (W{('-' if charge == 'minus' else '+')} {i}, {j})" for i in range(1, 6)
                       for j in range(1, 6) for charge in ["minus", "plus"]})
NP_DESCRIPTION.update({f"gamma_stat_{OSSS}_lep_{charge}_0tag_Dplus_{var}_bin{diff_bin}_bin_{mass_bin}": f"MC Stat. D+ ({OSSS} W{('-' if charge == 'minus' else '+')} {var}{diff_bin} mass{mass_bin})" for OSSS in [
    "OS", "SS"] for charge in ["minus", "plus"] for var in ["pt", "eta"] for diff_bin in range(1, 6) for mass_bin in range(0, 10)})
NP_DESCRIPTION.update({f"gamma_stat_{OSSS}_lep_{charge}_1tag_Dplus_bin_0": f"MC Stat. D+ ({OSSS} W{('-' if charge == 'minus' else '+')} 1-tag)" for OSSS in [
    "OS", "SS"] for charge in ["minus", "plus"]})
NP_DESCRIPTION.update({f"gamma_stat_{OSSS}_lep_{charge}_0tag_Dstar_{var}_bin{diff_bin}_bin_{mass_bin}": f"MC Stat. D* ({OSSS} W{('-' if charge == 'minus' else '+')} {var}{diff_bin} mass{mass_bin})" for OSSS in [
    "OS", "SS"] for charge in ["minus", "plus"] for var in ["pt", "eta"] for diff_bin in range(1, 6) for mass_bin in range(0, 13)})
NP_DESCRIPTION.update({f"gamma_stat_{OSSS}_lep_{charge}_1tag_Dstar_bin_0": f"MC Stat. D* ({OSSS} W{('-' if charge == 'minus' else '+')} 1-tag)" for OSSS in [
    "OS", "SS"] for charge in ["minus", "plus"]})

COMMON_NP_NAMES = {
    "EG_RESOLUTION_ALL": "EG_RESOLUTION_ALL",
    "EG_SCALE_ALL": "EG_SCALE_ALL",
    "EL_CHARGEID_STAT": "EL_CHARGEID_STAT",
    "EL_CHARGEID_SYStotal": "EL_CHARGEID_SYStotal",
    "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR": "EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR",
    "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR": "EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR",
    "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR": "EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR",
    "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR": "EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR",
    "FID_EFF_AS": "FID_EFF_AS",
    "FID_EFF_PDF": "FID_EFF_PDF",
    "FID_EFF_QCD": "FID_EFF_QCD",
    "FID_EFF_WPLUSD_MG": "FID_EFF_WPLUSD_MG",
    "JET_EffectiveNP_8restTerm": "JET_EffectiveNP_8restTerm",
    "JET_EtaIntercalibration_Modelling": "JET_EtaIntercalibration_Modelling",
    "JET_EtaIntercalibration_NonClosure_highE": "JET_EtaIntercalibration_NonClosure_highE",
    "JET_EtaIntercalibration_NonClosure_negEta": "JET_EtaIntercalibration_NonClosure_negEta",
    "JET_EtaIntercalibration_NonClosure_posEta": "JET_EtaIntercalibration_NonClosure_posEta",
    "JET_EtaIntercalibration_TotalStat": "JET_EtaIntercalibration_TotalStat",
    "JET_Flavor_Composition": "JET_Flavor_Composition",
    "JET_Flavor_Response": "JET_Flavor_Response",
    "JET_JER_DataVsMC_MC16": "JET_JER_DataVsMC_MC16",
    "JET_JER_EffectiveNP_7restTerm": "JET_JER_EffectiveNP_7restTerm",
    "JET_JvtEfficiency": "JET_JvtEfficiency",
    "JET_Pileup_OffsetMu": "JET_Pileup_OffsetMu",
    "JET_Pileup_OffsetNPV": "JET_Pileup_OffsetNPV",
    "JET_Pileup_PtTerm": "JET_Pileup_PtTerm",
    "JET_Pileup_RhoTopology": "JET_Pileup_RhoTopology",
    "JET_PunchThrough_MC16": "JET_PunchThrough_MC16",
    "Lumi": "Lumi",
    "MET_SoftTrk_ResoPara": "MET_SoftTrk_ResoPara",
    "MET_SoftTrk_ResoPerp": "MET_SoftTrk_ResoPerp",
    "MET_SoftTrk": "MET_SoftTrk",
    "mu_Top": "mu_Top",
    "MUON_CB": "MUON_CB",
    "MUON_EFF_ISO_STAT": "MUON_EFF_ISO_STAT",
    "MUON_EFF_ISO_SYS": "MUON_EFF_ISO_SYS",
    "MUON_EFF_RECO_STAT_LOWPT": "MUON_EFF_RECO_STAT_LOWPT",
    "MUON_EFF_RECO_STAT": "MUON_EFF_RECO_STAT",
    "MUON_EFF_RECO_SYS_LOWPT": "MUON_EFF_RECO_SYS_LOWPT",
    "MUON_EFF_RECO_SYS": "MUON_EFF_RECO_SYS",
    "MUON_EFF_TrigStatUncertainty": "MUON_EFF_TrigStatUncertainty",
    "MUON_EFF_TrigSystUncertainty": "MUON_EFF_TrigSystUncertainty",
    "MUON_EFF_TTVA_STAT": "MUON_EFF_TTVA_STAT",
    "MUON_EFF_TTVA_SYS": "MUON_EFF_TTVA_SYS",
    "MUON_SAGITTA_DATASTAT": "MUON_SAGITTA_DATASTAT",
    "MUON_SAGITTA_RESBIAS": "MUON_SAGITTA_RESBIAS",
    "MUON_SCALE": "MUON_SCALE",
    "PROD_FRAC_EIG_1": "PROD_FRAC_EIG_1",
    "PROD_FRAC_EIG_2": "PROD_FRAC_EIG_2",
    "PROD_FRAC_EIG_3": "PROD_FRAC_EIG_3",
    "PRW_DATASF": "PRW_DATASF",
    "Top_GEN_isr": "Top_GEN_isr",
    "Top_GEN_muF": "Top_GEN_muF",
    "Top_GEN_muR": "Top_GEN_muR",
    "Top_GEN_PDF": "Top_GEN_PDF",
    "Top_GEN_Var3c": "Top_GEN_Var3c",
    "Top_HS": "Top_HS",
    "Top_Shower": "Top_Shower",
    "TRK_EFF_D0_DEAD": "TRK_EFF_D0_DEAD",
    "TRK_EFF_D0_MEAS": "TRK_EFF_D0_MEAS",
    "TRK_EFF_IBL": "TRK_EFF_IBL",
    "TRK_EFF_Overal": "TRK_EFF_Overal",
    "TRK_EFF_PP0": "TRK_EFF_PP0",
    "TRK_EFF_QGSP": "TRK_EFF_QGSP",
    "TRK_EFF_Z0_DEAD": "TRK_EFF_Z0_DEAD",
    "TRK_EFF_Z0_MEAS": "TRK_EFF_Z0_MEAS",
    "TRK_EFF_D0_DEAD_Top": "TRK_EFF_D0_DEAD_Top",
    "TRK_EFF_D0_MEAS_Top": "TRK_EFF_D0_MEAS_Top",
    "TRK_EFF_IBL_Top": "TRK_EFF_IBL_Top",
    "TRK_EFF_Overal_Top": "TRK_EFF_Overal_Top",
    "TRK_EFF_PP0_Top": "TRK_EFF_PP0_Top",
    "TRK_EFF_QGSP_Top": "TRK_EFF_QGSP_Top",
    "TRK_EFF_Z0_DEAD_Top": "TRK_EFF_Z0_DEAD_Top",
    "TRK_EFF_Z0_MEAS_Top": "TRK_EFF_Z0_MEAS_Top",
}
COMMON_NP_NAMES.update({f"FT_EFF_Eigen_B_{i}": f"FT_EFF_Eigen_B_{i}" for i in range(0, 9)})
COMMON_NP_NAMES.update({f"FT_EFF_Eigen_C_{i}": f"FT_EFF_Eigen_C_{i}" for i in range(0, 9)})
COMMON_NP_NAMES.update({f"FT_EFF_Eigen_Light_{i}": f"FT_EFF_Eigen_Light_{i}" for i in range(0, 9)})
COMMON_NP_NAMES.update({f"JET_EffectiveNP_{i}": f"JET_EffectiveNP_{i}" for i in range(1, 10)})
COMMON_NP_NAMES.update({f"JET_JER_EffectiveNP_{i}": f"JET_JER_EffectiveNP_{i}" for i in range(1, 7)})

NP_NAMES = {
    "dplus": deepcopy(COMMON_NP_NAMES),
    "dstar": deepcopy(COMMON_NP_NAMES),
}

# D+ only uncertainites
NP_NAMES["dplus"].update({
    "BranchingRatio": "DPLUS_BR_SIG",
    "DMESON_MASS_TOT_minus": "DPLUS_PEAK_POSITION_PLUS",
    "DMESON_MASS_TOT_plus": "DPLUS_PEAK_POSITION_MINUS",
    "DMESON_RESO_TOT_minus": "DPLUS_PEAK_RESOLUTION_PLUS",
    "DMESON_RESO_TOT_plus": "DPLUS_PEAK_RESOLUTION_MINUS",
    "DPLUS_BKG_BR": "DPLUS_BR_BKG",
    "MM_EL_FR__Wjets_MadGraph": "DPLUS_MM_EL_FR_Wjets_MadGraph",
    "MM_EL_FR_MET_Var": "DPLUS_MM_EL_FR_MET_Var",
    "MM_EL_FR_STAT": "DPLUS_MM_EL_FR_STAT",
    "MM_EL_RR_MadGraph_Var": "DPLUS_MM_EL_RR_MadGraph_Var",
    "MM_EL_RR_STAT": "DPLUS_MM_EL_RR_STAT",
    "MM_MU_FR__Wjets_MadGraph": "DPLUS_MM_MU_FR_Wjets_MadGraph",
    "MM_MU_FR_MET_Var": "DPLUS_MM_MU_FR_MET_Var",
    "MM_MU_FR_STAT": "DPLUS_MM_MU_FR_STAT",
    "MM_MU_RR_MadGraph_Var": "DPLUS_MM_MU_RR_MadGraph_Var",
    "MM_MU_RR_STAT": "DPLUS_MM_MU_RR_STAT",
    "SHERPA2211_AS_Wcharm": "DPLUS_SHERPA2211_AS_WCHARM",
    "SHERPA2211_AS_Wjets": "DPLUS_SHERPA2211_AS_WJETS",
    "SHERPA2211_PDF_Wcharm": "DPLUS_SHERPA2211_PDF_WCHARM",
    "SHERPA2211_PDF_Wjets": "DPLUS_SHERPA2211_PDF_WJETS",
    "SHERPA2211_QCD_Wcharm": "DPLUS_SHERPA2211_QCD_WCHARM",
    "SHERPA2211_QCD_Wjets": "DPLUS_SHERPA2211_QCD_WJETS",
    "Matched_Sherpa_MG_2P": "DPLUS_Matched_Sherpa_MG_2P",
})

# D* only uncertainites
NP_NAMES["dstar"].update({
    "BranchingRatio": "DSTAR_BR_SIG",
    "DMESON_MASS_ADHOC_minus": "DSTAR_PEAK_POSITION_ADHOC_PLUS",
    "DMESON_MASS_ADHOC_plus": "DSTAR_PEAK_POSITION_ADHOC_MINUS",
    "DMESON_MASS_TOT_minus": "DSTAR_PEAK_POSITION_PLUS",
    "DMESON_MASS_TOT_plus": "DSTAR_PEAK_POSITION_MINUS",
    "DMESON_RESO_TOT_minus": "DSTAR_PEAK_RESOLUTION_PLUS",
    "DMESON_RESO_TOT_plus": "DSTAR_PEAK_RESOLUTION_MINUS",
    "MM_EL_FR__Wjets_MadGraph": "DSTAR_MM_EL_FR_Wjets_MadGraph",
    "MM_EL_FR_MET_Var": "DSTAR_MM_EL_FR_MET_Var",
    "MM_EL_FR_STAT": "DSTAR_MM_EL_FR_STAT",
    "MM_EL_RR_MadGraph_Var": "DSTAR_MM_EL_RR_MadGraph_Var",
    "MM_EL_RR_STAT": "DSTAR_MM_EL_RR_STAT",
    "MM_MU_FR__Wjets_MadGraph": "DSTAR_MM_MU_FR_Wjets_MadGraph",
    "MM_MU_FR_MET_Var": "DSTAR_MM_MU_FR_MET_Var",
    "MM_MU_FR_STAT": "DSTAR_MM_MU_FR_STAT",
    "MM_MU_RR_MadGraph_Var": "DSTAR_MM_MU_RR_MadGraph_Var",
    "MM_MU_RR_STAT": "DSTAR_MM_MU_RR_STAT",
    "PRW_DATASF_WJETS": "PRW_DATASF_WJETS",
    "SHERPA2211_AS_WJETS": "DSTAR_SHERPA2211_AS_WJETS",
    "SHERPA2211_AS_WPLUSC": "DSTAR_SHERPA2211_AS_WPLUSC",
    "SHERPA2211_PDF_WJETS": "DSTAR_SHERPA2211_PDF_WJETS",
    "SHERPA2211_PDF_WPLUSC": "DSTAR_SHERPA2211_PDF_WPLUSC",
    "SHERPA2211_QCD_WJETS": "DSTAR_SHERPA2211_QCD_WJETS",
    "SHERPA2211_QCD_WPLUSC": "DSTAR_SHERPA2211_QCD_WPLUSC",
    "WJETS_FIT": "WJETS_FIT",
})

# bkg shape and normalization
# D+
NP_NAMES["dplus"].update({"CharmNorm_1tag": "DPLUS_CharmNorm_1tag"})
NP_NAMES["dplus"].update({"MisMatchNorm_1tag": "DPLUS_MisMatchNorm_1tag"})
NP_NAMES["dplus"].update({"OtherNorm_0tag": "DPLUS_OtherNorm_0tag"})
NP_NAMES["dplus"].update({"OtherNorm_1tag": "DPLUS_OtherNorm_1tag"})
NP_NAMES["dplus"].update({"WjetsNorm_1tag": "DPLUS_WjetsNorm_1tag"})
NP_NAMES["dplus"].update({f"CharmNorm_0tag_bin{i}": f"DPLUS_CharmNorm_0tag_bin{i}" for i in range(1, 6)})
NP_NAMES["dplus"].update({f"MisMatch_Sh_MG_2P_bin{i}": f"DPLUS_MisMatch_Sh_MG_2P_bin{i}" for i in range(1, 6)})
NP_NAMES["dplus"].update({f"MisMatchNorm_0tag_bin{i}": f"DPLUS_MisMatchNorm_0tag_bin{i}" for i in range(1, 6)})
NP_NAMES["dplus"].update({f"Wjets_Sh_MG_2P_bin{i}": f"DPLUS_Wjets_Sh_MG_2P_bin{i}" for i in range(1, 6)})
NP_NAMES["dplus"].update({f"WjetsNorm_0tag_bin{i}": f"DPLUS_WjetsNorm_0tag_bin{i}" for i in range(1, 6)})

# D*
NP_NAMES["dstar"].update({"CharmNorm_1tag": "DSTAR_CharmNorm_1tag"})
NP_NAMES["dstar"].update({"MisMatchNorm_1tag": "DSTAR_MisMatchNorm_1tag"})
NP_NAMES["dstar"].update({"OtherNorm_0tag": "DSTAR_OtherNorm_0tag"})
NP_NAMES["dstar"].update({"OtherNorm_1tag": "DSTAR_OtherNorm_1tag"})
NP_NAMES["dstar"].update({"WjetsNorm_1tag": "DSTAR_WjetsNorm_1tag"})
NP_NAMES["dstar"].update({f"CharmNorm_0tag_bin{i}": f"DSTAR_CharmNorm_0tag_bin{i}" for i in range(1, 6)})
NP_NAMES["dstar"].update({f"MisMatched_Sherpa_MG_2P_bin{i}": f"DSTAR_MisMatched_Sherpa_MG_2P_bin{i}" for i in range(1, 6)})
NP_NAMES["dstar"].update({f"MisMatchNorm_0tag_bin{i}": f"DSTAR_MisMatchNorm_0tag_bin{i}" for i in range(1, 6)})
NP_NAMES["dstar"].update({f"Wjets_MG_Sherpa_2P_bin{i}": f"DSTAR_Wjets_MG_Sherpa_2P_bin{i}" for i in range(1, 6)})
NP_NAMES["dstar"].update({f"WjetsNorm_0tag_bin{i}": f"DSTAR_WjetsNorm_0tag_bin{i}" for i in range(1, 6)})

# bkg normalization stat. unc.
# D+
NP_NAMES["dplus"].update({f"Stat_{charge}_DibosonZjets_0tag_{OSSS}_{var}_bin{diff_bin}": f"DPLUS_Stat_{charge}_DibosonZjets_0tag_{OSSS}_{var}_bin{diff_bin}" for charge in [
                         "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_NAMES["dplus"].update({f"Stat_{charge}_WCharm_0tag_{OSSS}_{var}_bin{diff_bin}": f"DPLUS_Stat_{charge}_WCharm_0tag_{OSSS}_{var}_bin{diff_bin}" for charge in [
                         "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_NAMES["dplus"].update({f"Stat_{charge}_Wjets_0tag_{OSSS}_{var}_bin{diff_bin}": f"DPLUS_Stat_{charge}_Wjets_0tag_{OSSS}_{var}_bin{diff_bin}" for charge in [
                         "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_NAMES["dplus"].update({f"Stat_{charge}_WMisMatched_0tag_{OSSS}_{var}_bin{diff_bin}": f"DPLUS_Stat_{charge}_WMisMatched_0tag_{OSSS}_{var}_bin{diff_bin}" for charge in [
                         "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_NAMES["dplus"].update({f"Stat_{charge}_Top_0tag_{OSSS}_{var}_bin{diff_bin}": f"DPLUS_Stat_{charge}_Top_0tag_{OSSS}_{var}_bin{diff_bin}" for charge in [
                         "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})

# D*
NP_NAMES["dstar"].update({f"Stat_{charge}_DibosonZjets_0tag_{OSSS}_{var}_bin{diff_bin}": f"DSTAR_Stat_{charge}_DibosonZjets_0tag_{OSSS}_{var}_bin{diff_bin}" for charge in [
                         "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_NAMES["dstar"].update({f"Stat_{charge}_WCharm_0tag_{OSSS}_{var}_bin{diff_bin}": f"DSTAR_Stat_{charge}_WCharm_0tag_{OSSS}_{var}_bin{diff_bin}" for charge in [
                         "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_NAMES["dstar"].update({f"Stat_{charge}_Wjets_0tag_{OSSS}_{var}_bin{diff_bin}": f"DSTAR_Stat_{charge}_Wjets_0tag_{OSSS}_{var}_bin{diff_bin}" for charge in [
                         "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_NAMES["dstar"].update({f"Stat_{charge}_WMisMatched_0tag_{OSSS}_{var}_bin{diff_bin}": f"DSTAR_Stat_{charge}_WMisMatched_0tag_{OSSS}_{var}_bin{diff_bin}" for charge in [
                         "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})
NP_NAMES["dstar"].update({f"Stat_{charge}_Top_0tag_{OSSS}_{var}_bin{diff_bin}": f"DSTAR_Stat_{charge}_Top_0tag_{OSSS}_{var}_bin{diff_bin}" for charge in [
                         "minus", "plus"] for OSSS in ["OS", "SS"] for var in ["pt", "eta"] for diff_bin in range(1, 6)})

# fiducial efficiency stat. unc.
NP_NAMES["dplus"].update({f"MC_Stat_Fid_Eff_{charge}_{i}_{j}": f"DPLUS_MC_Stat_Fid_Eff_{charge}_{i}_{j}" for i in range(1, 6)
                         for j in range(1, 6) for charge in ["minus", "plus"]})
NP_NAMES["dstar"].update({f"MC_Stat_Fid_Eff_{charge}_{i}_{j}": f"DSTAR_MC_Stat_Fid_Eff_{charge}_{i}_{j}" for i in range(1, 6)
                         for j in range(1, 6) for charge in ["minus", "plus"]})
# NP_NAMES["dplus"].update({f"MC_Stat_Fid_Eff_{charge}_{i}_{j}": f"Response Matrix Stat. D+ (W{('-' if charge == 'minus' else '+')} {i}, {j})" for i in range(1, 6)
#                          for j in range(1, 6) for charge in ["minus", "plus"]})
# NP_NAMES["dstar"].update({f"MC_Stat_Fid_Eff_{charge}_{i}_{j}": f"Response Matrix Stat. D* (W{('-' if charge == 'minus' else '+')} {i}, {j})" for i in range(1, 6)
#                          for j in range(1, 6) for charge in ["minus", "plus"]})

# gammas
# D+
NP_NAMES["dplus"].update({f"gamma_stat_{OSSS}_lep_{charge}_0tag_Dplus_{var}_bin{diff_bin}_bin_{mass_bin}": f"gamma_stat_{OSSS}_lep_{charge}_0tag_Dplus_{var}_bin{diff_bin}_bin_{mass_bin}" for OSSS in [
                         "OS", "SS"] for charge in ["minus", "plus"] for var in ["pt", "eta"] for diff_bin in range(1, 6) for mass_bin in range(0, 10)})
NP_NAMES["dplus"].update({f"gamma_stat_{OSSS}_lep_{charge}_1tag_Dplus_bin_0": f"gamma_stat_{OSSS}_lep_{charge}_1tag_Dplus_bin_0" for OSSS in [
                         "OS", "SS"] for charge in ["minus", "plus"]})
# NP_NAMES["dplus"].update({f"gamma_stat_{OSSS}_lep_{charge}_0tag_Dplus_{var}_bin{diff_bin}_bin_{mass_bin}": f"MC Stat. D+ ({OSSS} W{('-' if charge == 'minus' else '+')} {var}{diff_bin} mass{mass_bin})" for OSSS in [
#                          "OS", "SS"] for charge in ["minus", "plus"] for var in ["pt", "eta"] for diff_bin in range(1, 6) for mass_bin in range(0, 10)})
# NP_NAMES["dplus"].update({f"gamma_stat_{OSSS}_lep_{charge}_1tag_Dplus_bin_0": f"MC Stat. D+ ({OSSS} W{('-' if charge == 'minus' else '+')} 1-tag)" for OSSS in [
#                          "OS", "SS"] for charge in ["minus", "plus"]})

# D*
NP_NAMES["dstar"].update({f"gamma_stat_{OSSS}_lep_{charge}_0tag_Dstar_{var}_bin{diff_bin}_bin_{mass_bin}": f"gamma_stat_{OSSS}_lep_{charge}_0tag_Dstar_{var}_bin{diff_bin}_bin_{mass_bin}" for OSSS in [
                         "OS", "SS"] for charge in ["minus", "plus"] for var in ["pt", "eta"] for diff_bin in range(1, 6) for mass_bin in range(0, 13)})
NP_NAMES["dstar"].update({f"gamma_stat_{OSSS}_lep_{charge}_1tag_Dstar_bin_0": f"gamma_stat_{OSSS}_lep_{charge}_1tag_Dstar_bin_0" for OSSS in [
                         "OS", "SS"] for charge in ["minus", "plus"]})
# NP_NAMES["dstar"].update({f"gamma_stat_{OSSS}_lep_{charge}_0tag_Dstar_{var}_bin{diff_bin}_bin_{mass_bin}": f"MC Stat. D* ({OSSS} W{('-' if charge == 'minus' else '+')} {var}{diff_bin} mass{mass_bin})" for OSSS in [
#                          "OS", "SS"] for charge in ["minus", "plus"] for var in ["pt", "eta"] for diff_bin in range(1, 6) for mass_bin in range(0, 13)})
# NP_NAMES["dstar"].update({f"gamma_stat_{OSSS}_lep_{charge}_1tag_Dstar_bin_0": f"MC Stat. D* ({OSSS} W{('-' if charge == 'minus' else '+')} 1-tag)" for OSSS in [
#                          "OS", "SS"] for charge in ["minus", "plus"]})

# make output folder
if not os.path.isdir("hepdata"):
    os.makedirs("hepdata")


def atlas_rounding(central, up, dn=None, precision=None):
    if precision:
        central = round(central, precision)
        up = round(up, precision)
        if dn:
            dn = round(dn, precision)
        return central, up, dn
    up = round(up, -int(floor(log10(abs(up / 10.)))))
    if dn:
        dn = round(dn, -int(floor(log10(abs(dn / 10.)))))
        central = round(central, min(-int(floor(log10(abs(up / 10.)))), -int(floor(log10(abs(dn / 10.))))))
        return central, up, dn, min(-int(floor(log10(abs(up / 10.)))), -int(floor(log10(abs(dn / 10.)))))
    else:
        central = round(central, -int(floor(log10(abs(up / 10.)))))
        return central, up, -int(floor(log10(abs(up / 10.))))


def bin_edges(POI, var):
    if var == "eta":
        i = int(POI[-1])
        return (0.5 * (i - 1), 0.5 * i)


def dependable_dict(name):
    return {'header': {'name': name},
            'qualifiers': QUALIFIERS,
            'values': []}


def main():

    # main yaml file
    submission = [
        {
            "data_file": "dplus_minus.yaml",
            "description": "The 'OS-SS' W+D fiducial phase-space absolute differential cross-sections in the W-D+ channel.",
            "keywords": [
                {"name": "reactions", "values": ["P P --> W D"]},
                {"name": "observables", "values": ["SIG"]},
                {"name": "cmenergies", "values": [13000.0]},
                {"name": "phrases", "values": ["Proton-Proton Scattering", "W Production", "D production"]},
            ],
            "location": "XXX",
            "name": "Table 7"
        },
        {
            "data_file": "dplus_plus.yaml",
            "description": "The 'OS-SS' W+D fiducial phase-space absolute differential cross-sections in the W+D- channel.",
            "keywords": [
                {"name": "reactions", "values": ["P P --> W D"]},
                {"name": "observables", "values": ["SIG"]},
                {"name": "cmenergies", "values": [13000.0]},
                {"name": "phrases", "values": ["Proton-Proton Scattering", "W Production", "D production"]},
            ],
            "location": "XXX",
            "name": "Table 8"
        },
        {
            "data_file": "dstar_minus.yaml",
            "description": "The 'OS-SS' W+D* fiducial phase-space absolute differential cross-sections in the W-D*+ channel.",
            "keywords": [
                {"name": "reactions", "values": ["P P --> W D"]},
                {"name": "observables", "values": ["SIG"]},
                {"name": "cmenergies", "values": [13000.0]},
                {"name": "phrases", "values": ["Proton-Proton Scattering", "W Production", "D production"]},
            ],
            "location": "XXX",
            "name": "Table 9"
        },
        {
            "data_file": "dstar_plus.yaml",
            "description": "The 'OS-SS' W+D* fiducial phase-space absolute differential cross-sections in the W+D*- channel.",
            "keywords": [
                {"name": "reactions", "values": ["P P --> W D"]},
                {"name": "observables", "values": ["SIG"]},
                {"name": "cmenergies", "values": [13000.0]},
                {"name": "phrases", "values": ["Proton-Proton Scattering", "W Production", "D production"]},
            ],
            "location": "XXX",
            "name": "Table 10"
        },
    ]
    with open('hepdata/submission_temp.yaml', 'w') as yaml_file:
        for table in submission:
            yaml_file.write(f"# Nuisance parameter ranking for {table['data_file'].replace('.yaml', '')}\n")
            yaml_file.write("---\n")
            yaml.dump(table, yaml_file, default_flow_style=False)
            yaml_file.write("\n")

    # keep for all channels
    GLOBAL_NPs = []

    # list of all writen out NPs
    ALL_NPs = set()

    # loop first to extract the NPs
    RANKINGS = {}
    for channel in CHANNELS:

        RANKINGS[channel] = {}

        # folder
        if channel == "dplus":
            br = 1.0 * 2.0
            FOLDER = DPLUS_FOLDER
            meson_str = "Dplus"
        elif channel == "dstar":
            br = 0.677 * 2.0
            FOLDER = DSTAR_FOLDER
            meson_str = "Dstar"

        # fit results
        obs_fit_abs = f"WCharm_{meson_str}_lep_obs_OSSS_complete_xsec_alt_eta"
        f_result_abs = ROOT.TFile(os.path.join(FOLDER, obs_fit_abs, "Fits", f"{obs_fit_abs}.root"), "READ")
        f_result_abs_stat = ROOT.TFile(os.path.join(FOLDER, obs_fit_abs, "Fits", f"{obs_fit_abs}_statOnly.root"), "READ")
        fr_abs = f_result_abs.Get("nll_simPdf_newasimovData_with_constr")
        fr_abs_stat = f_result_abs_stat.Get("nll_simPdf_newasimovData_with_constr")

        # print
        for POI in POIs_abs:
            par = fr_abs.floatParsFinal().find(POI)
            par_stat = fr_abs_stat.floatParsFinal().find(POI)
            print(POI)
            print(par.getVal() / br, par.getErrorHi() / br, par.getErrorLo() / br)
            print(par_stat.getVal() / br, par_stat.getErrorHi() / br, par_stat.getErrorLo() / br)
            print(par.getVal() / br, (par.getErrorHi()**2 - par_stat.getErrorHi()**2)**0.5 / br, (par.getErrorLo()**2 - par_stat.getErrorLo()**2)**0.5 / br)

        # read rankings
        # dict_keys(['Name', 'NPhat', 'NPerrHi', 'NPerrLo', 'POIup', 'POIdown', 'POIupPreFit', 'POIdownPreFit'])
        NPs = []
        for POI in POIs_abs:

            with open(os.path.join(FOLDER, obs_fit_abs, f"Ranking_{POI}.yaml"), 'r') as stream:
                ranking = yaml.safe_load(stream)
                RANKINGS[channel][POI] = ranking
                print(f"Read ranking for {POI} {channel}, size: {len(ranking)}")
                assert len(ranking) > 100, "Too small!"

                # calculate total error from ranking
                err_up_total = 0
                err_dn_total = 0
                for NP in ranking:
                    for key in ["POIup", "POIdown"]:
                        if NP[key] > 0:
                            err_up_total += NP[key]**2
                        else:
                            err_dn_total += NP[key]**2
                err_up_total = err_up_total**0.5
                err_dn_total = err_dn_total**0.5
                print("initial total error for: ", channel, POI, err_up_total / br, f"-{err_dn_total / br}")

                # keep NPs adding up to some threshold
                err_up = 0
                err_dn = 0
                threshold_reached = False
                for NP in ranking:
                    if not threshold_reached:
                        NPs += [NP["Name"]]
                        for key in ["POIup", "POIdown"]:
                            if NP[key] > 0:
                                err_up += NP[key]**2
                            else:
                                err_dn += NP[key]**2
                    if (err_up**0.5 / err_up_total) > PRUNING_THRESHOLD and (err_dn**0.5 / err_dn_total) > PRUNING_THRESHOLD:
                        threshold_reached = True

        # add to list of global NPs
        GLOBAL_NPs += NPs

    # loop again to make the tables
    for channel in CHANNELS:

        # folder
        if channel == "dplus":
            FOLDER = DPLUS_FOLDER
            br = 1.0 * 2.0
            meson_str = "Dplus"
        elif channel == "dstar":
            FOLDER = DSTAR_FOLDER
            br = 0.677 * 2.0
            meson_str = "Dstar"

        # fit results
        obs_fit_abs = f"WCharm_{meson_str}_lep_obs_OSSS_complete_xsec_alt_eta"
        f_result_abs = ROOT.TFile(os.path.join(FOLDER, obs_fit_abs, "Fits", f"{obs_fit_abs}.root"), "READ")
        f_result_abs_stat = ROOT.TFile(os.path.join(FOLDER, obs_fit_abs, "Fits", f"{obs_fit_abs}_statOnly.root"), "READ")
        fr_abs = f_result_abs.Get("nll_simPdf_newasimovData_with_constr")
        fr_abs_stat = f_result_abs_stat.Get("nll_simPdf_newasimovData_with_constr")

        # list of NPs (same for all channels)
        NPs = list(dict.fromkeys(GLOBAL_NPs))
        assert "MUON_SAGITTA_RESBIAS" in NPs, "MUON_SAGITTA_RESBIAS not included.. increase the threshold!"

        # get the remaining error of the pruned away NPs
        for POI in POIs_abs:

            # load the ranking
            ranking = RANKINGS[channel][POI]

            # keep NPs adding up to some threshold
            err_up_remaining = 0
            err_dn_remaining = 0
            count = 0
            for NP in ranking:
                if NP['Name'] not in NPs:
                    count += 1
                    # print(f"Adding up remaining error for {NP['Name']} for {POI}")
                    for key in ["POIup", "POIdown"]:
                        if NP[key] > 0:
                            err_up_remaining += NP[key]**2
                        else:
                            err_dn_remaining += NP[key]**2
            print("partial error for: ", channel, POI, err_up_remaining**0.5, f"-{err_dn_remaining**0.5}", count)
            if abs(err_up_remaining) > 0 or abs(err_dn_remaining) > 0:
                RANKINGS[channel][POI] += [{
                    'Name': f"Uncorr_{channel}_{POI}",
                    'POIup': err_up_remaining**0.5,
                    'POIdown': -err_dn_remaining**0.5,
                }]

        # NP impact table
        for charge in ["minus", "plus"]:

            # POIs
            if charge == "plus":
                POIs_channel = POIs_abs[0:5]
            elif charge == "minus":
                POIs_channel = POIs_abs[5:10]

            impact_table = {
                'independent_variables': [
                    {'header': {'name': 'LEP_ABS_ETA'},
                     'values': [{'high': bin_edges(POI, "eta")[0], 'low': bin_edges(POI, "eta")[1]} for POI in POIs_channel]}
                ],
                'dependent_variables': [
                    dependable_dict("DSIG/DLEP_ABS_ETA"),
                    dependable_dict("DSIG/DLEP_ABS_ETA (breakdown of systematics)"),
                ]
            }

            # fill in POIs
            for i, POI in enumerate(POIs_channel):
                par = fr_abs.floatParsFinal().find(POI)
                par_stat = fr_abs_stat.floatParsFinal().find(POI)
                vals_stat = atlas_rounding(par_stat.getVal() / br, par_stat.getErrorHi() / br, par_stat.getErrorLo() / br)
                vals = atlas_rounding(par.getVal() / br, (par.getErrorHi()**2 - par_stat.getErrorHi()**2)**0.5 / br,
                                      (par.getErrorLo()**2 - par_stat.getErrorLo()**2)**0.5 / br, vals_stat[-1])

                # create the two columns (2nd one gets filled out later)
                impact_table['dependent_variables'][0]['values'] += [{
                    'errors': [{'label': 'stat', 'symerror': vals_stat[1]},
                               {'label': 'sys', 'asymerror': {'minus': -vals[2], 'plus': vals[1]}}],
                    'value': vals[0]
                }]
                impact_table['dependent_variables'][1]['values'] += [{
                    'errors': [{'label': 'stat', 'symerror': vals_stat[1]}],
                    'value': vals[0]
                }]

            # sort NPs for each channel
            NPs_sorted = []
            for NP in NPs:
                sum_error = 0
                for POI in POIs_channel:
                    ranking = RANKINGS[channel][POI]
                    NP_impact = [x for x in ranking if x['Name'] == NP]
                    if len(NP_impact):
                        sum_error += max(abs(float(NP_impact[0]['POIup'])), abs(float(NP_impact[0]['POIdown'])))
                NPs_sorted += [(sum_error, NP)]
            NPs_sorted.sort(key=lambda x: x[0], reverse=True)

            # add the total correlated errors
            for POI in POIs_channel:
                NPs_sorted += [(0.0, f"Uncorr_{channel}_{POI}")]

            # fill rankings
            # total error validation
            err_up = {POI: 0 for POI in POIs_channel}
            err_dn = {POI: 0 for POI in POIs_channel}
            for _, NP in NPs_sorted:

                # loop over POIs
                for i, POI in enumerate(POIs_channel):

                    # load the ranking yaml
                    ranking = RANKINGS[channel][POI]

                    # # skip mu_Top
                    # if NP == "mu_Top":
                    #     print("skipping")
                    #     continue

                    # extract NP impact
                    NP_impact = [x for x in ranking if x['Name'] == NP]
                    if len(NP_impact):
                        vals = atlas_rounding(par.getVal() / br, float(NP_impact[0]['POIup']) / br, float(NP_impact[0]['POIdown']) / br)
                    else:
                        # print(f"WARNING: NP {NP} not found for POI {POI} in channel {channel}")
                        continue

                    # add NP in table
                    if NP in NP_NAMES[channel]:
                        impact_table['dependent_variables'][1]['values'][i]['errors'] += [
                            {'label': NP_NAMES[channel][NP], 'asymerror': {'minus': vals[2], 'plus': vals[1]}}]
                        ALL_NPs.add(NP_NAMES[channel][NP])
                    else:
                        print("Missing map for NP: ", NP)
                        impact_table['dependent_variables'][1]['values'][i]['errors'] += [{'label': NP, 'asymerror': {'minus': vals[2], 'plus': vals[1]}}]
                        ALL_NPs.add(NP)
                    # impact_table['dependent_variables'][1]['values'][i]['errors'] += [{'label': NP, 'asymerror': {'minus': vals[2], 'plus': vals[1]}}]

                    # validation
                    for i in [1, 2]:
                        if vals[i] > 0:
                            err_up[POI] += vals[i]**2
                        else:
                            err_dn[POI] += vals[i]**2

            for POI in POIs_channel:
                print(f"final total error for {channel} {POI}: {err_up[POI]**0.5} -{err_dn[POI]**0.5}")

            with open(f'hepdata/{channel}_{charge}.yaml', 'w') as yaml_file:
                yaml.dump(impact_table, yaml_file, default_flow_style=False)

    for NP in sorted(ALL_NPs):
        print(f"  - {{value: '{NP}'}}")
    print("-----------------")
    for NP in sorted(ALL_NPs):
        if NP in NP_DESCRIPTION:
            print(f"  - {{value: '{NP_DESCRIPTION[NP]}'}}")
        else:
            print(f"  - {{value: '{NP}'}}")

    fout = open("np_names.txt", "w")
    for NP in sorted(ALL_NPs):
        NP_name = NP
        if NP_name.startswith("gamma_stat_"):
            NP_name = NP_name.replace("gamma_stat_", "$\\gamma$(")
            NP_name += ")"
        NP_name = NP_name.replace("_", "\\_")
        if "Uncorr_" in NP:
            meson = "D+" if "dplus" in NP else "D*"
            wboson = "W-" if "Wminus" in NP else "W+"
            string = f"{NP_name} & The rest of the error for {meson} ({wboson} bin {NP[-1]} and not correlated with other bins) \\\\ \\hline\n"
        else:
            string = f"{NP_name} & {NP_DESCRIPTION[NP]} \\\\ \\hline\n"
        fout.write(string)
    fout.close()

if __name__ == "__main__":

    # run
    main()
