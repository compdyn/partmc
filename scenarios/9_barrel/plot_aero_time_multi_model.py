#!/usr/bin/env python
import sys
sys.path.append('../../tool/')
import subprocess
import os
import math
import numpy as np
import mpl_helper
import matplotlib.pyplot as plt

barrel_time = np.loadtxt("ref_aerodyne_0326/ref_aero_time.txt")
partmc_time1 = np.loadtxt("out_aerodyne_0326/barrel_wc_0001_aero_time_kd_point005.txt")
partmc_time2 = np.loadtxt("out_aerodyne_0326/barrel_wc_0001_aero_time_kd_point01.txt")
partmc_time3 = np.loadtxt("out_aerodyne_0326/barrel_wc_0001_aero_time_kd_point05.txt")
(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.plot(barrel_time[:,0], barrel_time[:,1], color='k')
axes.plot(partmc_time1[:,0], partmc_time1[:,1], color='r')
axes.plot(partmc_time2[:,0], partmc_time2[:,1], color='b')
axes.plot(partmc_time3[:,0], partmc_time3[:,1], color='g')
axes.set_xlabel("Time (s)")
axes.set_ylabel(r"Number concentration ($\mathrm{m}^{-3}$)")
axes.grid()
#axes.set_yscale('log')
axes.legend(('Barrel', r'$k_{\rm D}$ = 0.005', r'$k_{\rm D}$ = 0.01', r'$k_{\rm D}$ = 0.05'))
filename_out = "aero_num_time.pdf"
figure.savefig(filename_out)
