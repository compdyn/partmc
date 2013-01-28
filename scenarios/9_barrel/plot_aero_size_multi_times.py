#!/usr/bin/env python
import sys
sys.path.append('../../tool/')
import subprocess
import os
import math
import numpy as np
import mpl_helper
import matplotlib.pyplot as plt

# plot number distribution
partmc_num = np.loadtxt("out_0322/barrel_wc_0001_aero_size_num.txt")
barrel_num = np.loadtxt("ref_0322/ref_aero_size_num_regrid.txt")
(figure, axes) = mpl_helper.make_fig(colorbar=False)
b1 = axes.semilogx(barrel_num[:,0], barrel_num[:,1], color='r')
axes.semilogx(partmc_num[:,0], partmc_num[:,1]*math.log(10), color='k',linestyle='--')
b2 = axes.semilogx(barrel_num[:,0], barrel_num[:,11], color='b')
axes.semilogx(partmc_num[:,0], partmc_num[:,11]*math.log(10), color='k',linestyle='--')
b3 = axes.semilogx(barrel_num[:,0], barrel_num[:,21], color='g')
axes.semilogx(partmc_num[:,0], partmc_num[:,21]*math.log(10), color='k',linestyle='--')
b4 = axes.semilogx(barrel_num[:,0], barrel_num[:,31], color='y')
axes.semilogx(partmc_num[:,0], partmc_num[:,31]*math.log(10), color='k',linestyle='--')
b5 = axes.semilogx(barrel_num[:,0], barrel_num[:,39], color='c')
p5 = axes.semilogx(partmc_num[:,0], partmc_num[:,39]*math.log(10), color='k',linestyle='--')
axes.set_xlabel("Diameter (m)")
axes.set_ylabel(r"Number concentration ($\mathrm{m}^{-3}$)")
axes.grid()
axes.set_ylim(0, 1.05*max(partmc_num[:,1]*math.log(10)))
axes.legend([b1,b2,b3,b4,b5,p5], ('0 min', '70 min', '140 min', \
                                  '210 min', '266 min', 'PartMC'), loc='upper left')
filename_out = "aero_num_size_multi_times.pdf"
figure.savefig(filename_out)

# plot mass distribution
rho = 1760 # density in kgm-3
(figure, axes) = mpl_helper.make_fig(colorbar=False)
b1 = axes.semilogx(barrel_num[:,0], barrel_num[:,1] * math.pi / 6. * rho * barrel_num[:,0]**3, color='r')
axes.semilogx(partmc_num[:,0], partmc_num[:,1] * math.pi / 6. * rho * partmc_num[:,0]**3 * math.log(10), color='k', linestyle='--')
b2 = axes.semilogx(barrel_num[:,0], barrel_num[:,11] * math.pi / 6. * rho * barrel_num[:,0]**3, color='b')
axes.semilogx(partmc_num[:,0], partmc_num[:,11] * math.pi / 6. * rho * partmc_num[:,0]**3 * math.log(10), color='k', linestyle='--')
b3 = axes.semilogx(barrel_num[:,0], barrel_num[:,21] * math.pi / 6. * rho * barrel_num[:,0]**3, color='g')
axes.semilogx(partmc_num[:,0], partmc_num[:,21] * math.pi / 6. * rho * partmc_num[:,0]**3 * math.log(10), color='k', linestyle='--')
b4 = axes.semilogx(barrel_num[:,0], barrel_num[:,31] * math.pi / 6. * rho * barrel_num[:,0]**3, color='y')
axes.semilogx(partmc_num[:,0], partmc_num[:,31] * math.pi / 6. * rho * partmc_num[:,0]**3 * math.log(10), color='k', linestyle='--')
b5 = axes.semilogx(barrel_num[:,0], barrel_num[:,39] * math.pi / 6. * rho * barrel_num[:,0]**3, color='c')
p5 = axes.semilogx(partmc_num[:,0], partmc_num[:,39] * math.pi / 6. * rho * partmc_num[:,0]**3 * math.log(10), color='k', linestyle='--')
axes.set_xlabel("Diameter (m)")
axes.set_ylabel(r"Mass concentration (kg $\mathrm{m}^{-3}$)")
axes.grid()
axes.set_ylim(0, 1.05*max(partmc_num[:,1] * math.pi / 6. * rho * partmc_num[:,0]**3 * math.log(10)))
axes.legend([b1,b2,b3,b4,b5,p5], ('0 min', '70 min', '140 min', \
                                  '210 min', '266 min', 'PartMC'), loc='upper left')
filename_out = "aero_mass_size_multi_times.pdf"
figure.savefig(filename_out)
