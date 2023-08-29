#!/usr/bin/env python

from visual_tools.case import PartMC_Case

a = PartMC_Case("sdm_a")
mixing, dim, info = a.compute_ice_ratio("0001")
print(dim)
print(info)
print(mixing)

import matplotlib.pyplot as plt
plt.plot(mixing)
plt.grid()
plt.show()
