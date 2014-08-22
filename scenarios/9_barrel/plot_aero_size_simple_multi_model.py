#!/usr/bin/env python
import sys
sys.path.append('../../tool/')
import subprocess
import os
import math
import numpy as np
import mpl_helper
import matplotlib.pyplot as plt

col = 121

partmc_num1 = np.loadtxt("out_aerodyne_0326/barrel_wc_0001_aero_size_num_kd_point005.txt")
partmc_num2 = np.loadtxt("out_aerodyne_0326/barrel_wc_0001_aero_size_num_kd_point01.txt")
partmc_num3 = np.loadtxt("out_aerodyne_0326/barrel_wc_0001_aero_size_num_kd_point05.txt")
barrel_num = np.loadtxt("ref_aerodyne_0326/ref_aero_size_num.txt")

(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.semilogx(barrel_num[:,0], barrel_num[:,col], color='k')
axes.semilogx(partmc_num1[:,0], partmc_num1[:,col]*math.log(10), color='r')
axes.semilogx(partmc_num2[:,0], partmc_num2[:,col]*math.log(10), color='b')
axes.semilogx(partmc_num3[:,0], partmc_num3[:,col]*math.log(10), color='g')
axes.set_xlabel("Diameter (m)")
axes.set_ylabel(r"Number concentration ($\mathrm{m}^{-3}$)")
axes.grid()
axes.set_ylim(0, 1.1*max(partmc_num3[:,col])*math.log(10))
axes.legend(('Barrel', r'$k_{\rm D}$ = 0.005', r'$k_{\rm D}$ = 0.01', r'$k_{\rm D}$ = 0.05'), loc='best')
#axes.legend(('Barrel', r'$d_{\rm f}$ = 2', r'$d_{\rm f}$ = 2.3', r'$d_{\rm f}$ = 3'), loc='best')
filename_out = "aero_num_size.pdf"
figure.savefig(filename_out)
