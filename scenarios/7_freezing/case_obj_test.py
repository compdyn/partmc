#!/usr/bin/env python

from visual_tools.case import PartMC_Case
import numpy as np
from pdb import set_trace

#a = PartMC_Case("out_sdm_a", OutDir = "output")
a = PartMC_Case("pap1_bin_chis_constT_0.5", OutDir = "/data/keeling/a/wenhant2/d/modeldata/partmc_cases/pap1_bin_chis_10", disable_check_time = True)
ice_ratio, dim, info = a.compute_ice_ratio("0001")
#print(dim)
#print(info)
#print(ice_ratio)
for ensem in ["0001", "0002", "0003"]:
    Dp = a.compute_dry_diameter("0001")
    print(10**np.mean(np.log10(Dp[0])))
    print(np.std(np.log10(Dp[0])))

#import matplotlib.pyplot as plt
#plt.plot(ice_ratio)
#plt.grid()
#plt.show()
