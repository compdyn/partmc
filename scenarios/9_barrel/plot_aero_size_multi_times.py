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
partmc_num = np.loadtxt("out_aerodyne_0716/barrel_wc_0001_aero_size_num.txt")
barrel_num = np.loadtxt("ref_aerodyne_0716/ref_aero_size_num.txt")

(figure, axes) = mpl_helper.make_fig(colorbar=False)
b1, = axes.semilogx(barrel_num[:,0], barrel_num[:,1], color='r')
axes.semilogx(partmc_num[:,0], partmc_num[:,1]*math.log(10), color='k')
b2, = axes.semilogx(barrel_num[:,0], barrel_num[:,21], color='b')
axes.semilogx(partmc_num[:,0], partmc_num[:,21]*math.log(10), color='k')
b3, = axes.semilogx(barrel_num[:,0], barrel_num[:,41], color='g')
axes.semilogx(partmc_num[:,0], partmc_num[:,41]*math.log(10), color='k')
b4, = axes.semilogx(barrel_num[:,0], barrel_num[:,81], color='y')
axes.semilogx(partmc_num[:,0], partmc_num[:,81]*math.log(10), color='k')
b5, = axes.semilogx(barrel_num[:,0], barrel_num[:,121], color='c')
p5, = axes.semilogx(partmc_num[:,0], partmc_num[:,121]*math.log(10), color='k')
axes.set_xlabel("Diameter (m)")
axes.set_ylabel(r"Number concentration ($\mathrm{m}^{-3}$)")
axes.grid()
axes.set_ylim(0, 1.05*max(partmc_num[:,1]*math.log(10)))
#axes.set_ylim(0,2.6e11)
axes.legend([b1,b2,b3,b4,b5,p5], ('0 min', '40 min', '80 min', \
                                  '160 min', '240 min', 'PartMC'), loc='best')
#bbox_props = dict(boxstyle="square,pad=0.3", fc="cyan", ec="b", lw=1)
#axes.annotate('after 1hr', xy=(0.05, 0.9), xycoords='axes fraction',weight='extra bold', size=14, bbox=bbox_props)
filename_out = "aero_num_size_multi_times.pdf"
figure.savefig(filename_out)

# plot mass distribution
#rho = 1760 # density in kgm-3
#(figure, axes) = mpl_helper.make_fig(colorbar=False)
#b1 = axes.semilogx(barrel_num[:,0], barrel_num[:,1] * math.pi / 6. * rho * barrel_num[:,0]**3, color='r')
#axes.semilogx(partmc_num[:,0], partmc_num[:,1] * math.pi / 6. * rho * partmc_num[:,0]**3 * math.log(10), color='k', linestyle='--')
#b2 = axes.semilogx(barrel_num[:,0], barrel_num[:,11] * math.pi / 6. * rho * barrel_num[:,0]**3, color='b')
#axes.semilogx(partmc_num[:,0], partmc_num[:,11] * math.pi / 6. * rho * partmc_num[:,0]**3 * math.log(10), color='k', linestyle='--')
#b3 = axes.semilogx(barrel_num[:,0], barrel_num[:,21] * math.pi / 6. * rho * barrel_num[:,0]**3, color='g')
#axes.semilogx(partmc_num[:,0], partmc_num[:,21] * math.pi / 6. * rho * partmc_num[:,0]**3 * math.log(10), color='k', linestyle='--')
#b4 = axes.semilogx(barrel_num[:,0], barrel_num[:,31] * math.pi / 6. * rho * barrel_num[:,0]**3, color='y')
#axes.semilogx(partmc_num[:,0], partmc_num[:,31] * math.pi / 6. * rho * partmc_num[:,0]**3 * math.log(10), color='k', linestyle='--')
#b5 = axes.semilogx(barrel_num[:,0], barrel_num[:,41] * math.pi / 6. * rho * barrel_num[:,0]**3, color='c')
#p5 = axes.semilogx(partmc_num[:,0], partmc_num[:,41] * math.pi / 6. * rho * partmc_num[:,0]**3 * math.log(10), color='k', linestyle='--')
#axes.set_xlabel("Diameter (m)")
#axes.set_ylabel(r"Mass concentration (kg $\mathrm{m}^{-3}$)")
#axes.grid()
#axes.set_ylim(0, 1.05*max(partmc_num[:,11] * math.pi / 6. * rho * partmc_num[:,0]**3 * math.log(10)))
#axes.set_ylim(0, 8e-7)
#axes.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
#axes.legend([b1,b2,b3,b4,b5,p5], ('0 min', '70 min', '140 min', \
#                                  '210 min', '280 min', 'PartMC'), loc='upper left')
#axes.annotate('$d_f$ = 2.4', xy=(0.8, 0.9), xycoords='axes fraction',weight='extra bold', size=14, bbox=bbox_props)
#filename_out = "aero_mass_size_multi_times.pdf"
#figure.savefig(filename_out)
