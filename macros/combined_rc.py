# D+
# 0.9650511142497966 +0.006742196904525775 -0.006695226662652633 (stat)
# 0.9650511142497966 +0.013464383653938029 -0.013305555331056327 (sys + stat)
# 0.9650511142497966 +0.011654716130435637 -0.011498336514622743 (sys only)

# D*
# 0.9800957389237532 +0.009923325159115676 -0.009790778930560165 (stat)
# 0.9800957389237532 +0.016702717058781975 -0.016460894195604436 (sys + stat)
# 0.9800957389237532 +0.013435340521631308 -0.013232599353557976 (sys only)

# muon sagitta bias is a 0.3% uncertainty
muon_sagitta = 0.3 / 100.

dplus_central = 0.9650511142497966
dplus_stat = (0.006742196904525775 + 0.006695226662652633) / 2.
dplus_sys = (0.011654716130435637 + 0.011498336514622743) / 2.
dplus_sys_no_muon = (dplus_sys**2 - (muon_sagitta * dplus_central)**2)**0.5
dplus_total = (dplus_sys_no_muon**2 + dplus_stat**2)**0.5

dstar_central = 0.9800957389237532
dstar_stat = (0.009923325159115676 + 0.009790778930560165) / 2.
dstar_sys = (0.013435340521631308 + 0.013232599353557976) / 2.
dstar_sys_no_muon = (dstar_sys**2 - (muon_sagitta * dstar_central)**2)**0.5
dstar_total = (dstar_sys_no_muon**2 + dstar_stat**2)**0.5

# Eq. 40.8 in https://pdg.lbl.gov/2022/reviews/rpp2022-rev-statistics.pdf
w1 = 1 / (dplus_total**2)
w2 = 1 / (dstar_total**2)
w = w1 + w2

combined_central = (w1 * dplus_central + w2 * dstar_central) / w
combined_stat = 1 / ((1 / dplus_stat**2) + (1 / dstar_total**2))**0.5
combined_sys = 1 / w**0.5
combined_sys = (combined_sys**2 + (muon_sagitta * combined_central)**2)**0.5


print("Combined Rc")
print(f"central: {combined_central}")
print(f"stat: {combined_stat}")
print(f"sys: {combined_sys}")

# Combined Rc
# central: 0.9709292420861776
# stat: 0.006212765099712398
# sys: 0.010608271899824651
