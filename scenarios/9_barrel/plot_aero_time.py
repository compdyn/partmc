#!/usr/bin/env python
import sys
sys.path.append('../../tool/')
import subprocess
import os
import math
import numpy as np
import mpl_helper
import matplotlib.pyplot as plt

partmc_time = np.loadtxt("out_with_fractal/case_0008_wc_0001_aero_time.txt")
barrel_time = np.loadtxt("ref_aero_time.txt")
(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.plot(barrel_time[:,0], barrel_time[:,1], color='r')
axes.plot(partmc_time[:,0], partmc_time[:,1], color='k')
axes.set_xlabel("Time (s)")
axes.set_ylabel(r"Number concentration ($\mathrm{m}^{-3}$)")
axes.grid()
axes.legend(('Barrel', 'PartMC'))
filename_out = "aero_num_time.pdf"
figure.savefig(filename_out)
(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.plot(barrel_time[:,0], barrel_time[:,2], color='r')
axes.plot(partmc_time[:,0], partmc_time[:,2], color='k')
axes.set_xlabel("Time (s)")
axes.set_ylabel(r"Mass concentration (kg $\mathrm{m}^{-3}$)")
axes.grid()
axes.legend(('Barrel', 'PartMC'))
filename_out = "aero_mass_time.pdf"
figure.savefig(filename_out)
