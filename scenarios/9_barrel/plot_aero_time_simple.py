#!/usr/bin/env python
import sys
sys.path.append('../../tool/')
import subprocess
import os
import math
import numpy as np
import mpl_helper
import matplotlib.pyplot as plt

partmc_time = np.loadtxt("out_aerodyne_0325/barrel_wc_0001_aero_time.txt")
barrel_time = np.loadtxt("ref_aerodyne_0325/ref_aero_time.txt")
(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.plot(barrel_time[:,0], barrel_time[:,1], color='r')
axes.plot(partmc_time[:,0], partmc_time[:,1], color='k')
axes.set_xlabel("Time (s)")
axes.set_ylabel(r"Number concentration ($\mathrm{m}^{-3}$)")
axes.grid()
#axes.set_yscale('log')
axes.legend(('Barrel', 'PartMC'))
filename_out = "aero_num_time.pdf"
figure.savefig(filename_out)
