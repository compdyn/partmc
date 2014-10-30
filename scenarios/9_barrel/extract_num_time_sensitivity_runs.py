#!/usr/bin/env python
import sys
sys.path.append('../../tool/')
import subprocess
import os
import math
import numpy as np
import mpl_helper
import matplotlib.pyplot as plt

partmc_only_coag = np.loadtxt("out_0925/only_coag_0001_aero_time.txt")
partmc_only_wall = np.loadtxt("out_0925/only_wall_0001_aero_time.txt")
partmc_coag_wall = np.loadtxt("out_0925/case_0415_wc_0001_aero_time.txt")
(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.plot(partmc_only_coag[:,0], partmc_only_coag[:,1], color='r')
axes.plot(partmc_only_wall[:,0], partmc_only_wall[:,1], color='b')
axes.plot(partmc_coag_wall[:,0], partmc_coag_wall[:,1], color='k')
axes.set_xlabel("Time (s)")
axes.set_ylabel(r"Number concentration ($\mathrm{m}^{-3}$)")
axes.grid()
#axes.set_yscale('log')
axes.legend(('Only coag.', 'Only wall', 'Coag.+wall'), loc='best')
filename_out = "aero_num_time.pdf"
figure.savefig(filename_out)
