#!/usr/bin/env python
import sys
sys.path.append('../../tool/')
import subprocess
import os
import math
import numpy as np
import mpl_helper
import matplotlib.pyplot as plt

col = 11

partmc_num = np.loadtxt("out_aerodyne_0325/barrel_wc_0001_aero_size_num.txt")
barrel_num = np.loadtxt("ref_aerodyne_0325/ref_aero_size_num.txt")

(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.semilogx(barrel_num[:,0], barrel_num[:,col], color='r')
axes.semilogx(partmc_num[:,0], partmc_num[:,col]*math.log(10), color='k')
axes.set_xlabel("Diameter (m)")
axes.set_ylabel(r"Number concentration ($\mathrm{m}^{-3}$)")
axes.grid()
axes.set_ylim(0, 1.05*max(barrel_num[:,1]))
axes.legend(('Barrel', 'PartMC'))
filename_out = "aero_num_size.pdf"
figure.savefig(filename_out)
