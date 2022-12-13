# D+
# 0.964739849325106 +0.006729660597977746 -0.006685480592483229 (stat)
# 0.964739849325106 +0.013346278857945246 -0.013191764098044033 (sys + stat)
# 0.964739849325106 +0.011525399237773157 -0.011372202480872982 (sys only)


# D*
# 0.9797502160453198 +0.009917807554477426 -0.009782700185892685 (stat)
# 0.9797502160453198 +0.01668711742718649 -0.01644454866109697 (sys + stat)
# 0.9797502160453198 +0.013420021659485557 -0.013218243368168152 (sys only)

# muon sagitta bias is a 0.3% uncertainty
muon_sagitta = 0.3 / 100.

dplus_central = 0.964739849325106
dplus_stat = (0.006729660597977746 + 0.006685480592483229) / 2.
dplus_sys = (0.011525399237773157 + 0.011372202480872982) / 2.
dplus_sys_no_muon = (dplus_sys**2 - (muon_sagitta * dplus_central)**2)**0.5
dplus_total = (dplus_sys_no_muon**2 + dplus_stat**2)**0.5

dstar_central = 0.9797502160453198
dstar_stat = (0.009917807554477426 + 0.009782700185892685) / 2.
dstar_sys = (0.013420021659485557 + 0.013218243368168152) / 2.
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
# central: 0.9706063555826351
# stat: 0.006217251487304815
# sys: 0.010356382050784051
