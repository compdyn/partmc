#!/usr/bin/env python
import sys
sys.path.append('../../tool/')
import subprocess
import os
import math
import numpy as np
import mpl_helper
import matplotlib.pyplot as plt

prefactor = float(raw_input("Enter prefactor:"))
exponent = float(raw_input("Enter exponent:"))
data = np.loadtxt("rmse_total.dat")
for row in range(0, data.shape[0]):
    if data[row,1] == prefactor and data[row,2] == exponent:
       caseID = int(data[row,0])

filename_partmc_time = "out/case_%04d_wc_0001_aero_time.txt" %(caseID)

# extract aero time
if not os.path.isfile(filename_partmc_time):
   str_extr = "../../build/extract_aero_time"
   ncfile_prefix = "out/case_%04d_wc_0001" % (caseID)
   command = [str_extr, ncfile_prefix]
   subprocess.check_call(command)

# plot time series of number and mass concentrations
partmc_time = np.loadtxt(filename_partmc_time)
barrel_time = np.loadtxt("ref_aero_time.txt")
(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.plot(barrel_time[:,0], barrel_time[:,1])
axes.plot(partmc_time[:,0], partmc_time[:,1], marker='o', mfc='None')
axes.set_title("prefactor = %.3f, exponent = %.3f" % (prefactor, exponent))
axes.set_xlabel("Time (s)")
axes.set_ylabel("Number concentration ($\mathrm{m}^{-3}$)")
axes.grid()
axes.legend(('Barrel', 'PartMC'))
filename_out = "plot_aero_time_num.pdf"
figure.savefig(filename_out)
(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.plot(barrel_time[:,0], barrel_time[:,2])
axes.plot(partmc_time[:,0], partmc_time[:,2], marker='o', mfc='None')
axes.set_title("prefactor = %.3f, exponent = %.3f" % (prefactor, exponent))
axes.set_xlabel("Time (s)")
axes.set_ylabel("Mass concentration (kg $\mathrm{m}^{-3}$)")
axes.grid()
axes.legend(('Barrel', 'PartMC'))
filename_out = "plot_aero_time_mass.pdf"
figure.savefig(filename_out)

# plot number and mass distributions
partmc_num = np.loadtxt("out/case_%04d_wc_0001_aero_size_num.txt" %(caseID))
barrel_num = np.loadtxt("ref_aero_size_num_regrid.txt")
partmc_mass = np.loadtxt("out/case_%04d_wc_0001_aero_size_mass.txt" %(caseID))
barrel_mass = np.loadtxt("ref_aero_size_mass_regrid.txt")
for time in [0,10,20,30]:
    (figure, axes) = mpl_helper.make_fig(colorbar=False)
    axes.semilogx(barrel_num[:,0], barrel_num[:,time+1])
    axes.semilogx(partmc_num[:,0], partmc_num[:,time+1]*math.log(10), marker='o', mfc='None')
    axes.set_title("k\_{D} = %.3f, a = %.3f at time %02d" \
    % (prefactor, exponent, time))
    axes.set_xlabel("Dry diameter (m)")
    axes.set_ylabel("Number concentration ($\mathrm{m}^{-3}$)")
    axes.grid()
    axes.legend(('Barrel', 'PartMC'))
    filename_out = "plot_aero_size_num_time%02d.pdf" %(time)
    figure.savefig(filename_out)
for time in [0,10,20,30]:
    (figure, axes) = mpl_helper.make_fig(colorbar=False)
    axes.semilogx(barrel_mass[:,0], barrel_mass[:,time+1])
    axes.semilogx(partmc_mass[:,0], partmc_mass[:,time+1]*math.log(10), marker='o', mfc='None')
    axes.set_title("k\_{D} = %.3f, a = %.3f at time %02d" \
    % (prefactor, exponent, time))
    axes.set_xlabel("Dry diameter (m)")
    axes.set_ylabel("Mass concentration (kg $\mathrm{m}^{-3}$)")
    axes.grid()
    axes.legend(('Barrel', 'PartMC'), loc='upper left')
    filename_out = "plot_aero_size_mass_time%02d.pdf" %(time)
    figure.savefig(filename_out)
