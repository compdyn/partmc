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

max_time_num = []
min_time_num = []
max_time_mass = []
min_time_mass = []
time = []

data1 = np.loadtxt("out_with_bounds/0925_low_size_high_conc_0001_aero_time.txt")
data2 = np.loadtxt("out_with_bounds/0925_low_size_low_conc_0001_aero_time.txt")
data3 = np.loadtxt("out_with_bounds/0925_high_size_high_conc_0001_aero_time.txt")
data4 = np.loadtxt("out_with_bounds/0925_high_size_low_conc_0001_aero_time.txt")

for row in range(0,data1.shape[0]):
    max_time_num.append(max(data1[row,1],data2[row,1],data3[row,1],data4[row,1]))
    min_time_num.append(min(data1[row,1],data2[row,1],data3[row,1],data4[row,1]))
    max_time_mass.append(max(data1[row,2],data2[row,2],data3[row,2],data4[row,2]))
    min_time_mass.append(min(data1[row,2],data2[row,2],data3[row,2],data4[row,2]))
    time.append(data1[row,0])

array_max_size_num = np.loadtxt("out_with_bounds/0925_low_size_high_conc_0001_aero_size_num.txt")
array_min_size_num = np.loadtxt("out_with_bounds/0925_low_size_high_conc_0001_aero_size_num.txt")
array_max_size_mass = np.loadtxt("out_with_bounds/0925_low_size_high_conc_0001_aero_size_mass.txt")
array_min_size_mass = np.loadtxt("out_with_bounds/0925_low_size_high_conc_0001_aero_size_mass.txt")

for prefix in ['0925_low_size_low_conc', '0925_high_size_high_conc','0925_high_size_low_conc']:
    size_num = np.loadtxt("out_with_bounds/"+prefix+"_0001_aero_size_num.txt")
    for row in range(0,size_num.shape[0]):
        for col in range(1,size_num.shape[1]):
            array_max_size_num[row,col] = max(array_max_size_num[row,col],size_num[row,col])
            array_min_size_num[row,col] = min(array_min_size_num[row,col],size_num[row,col])
    size_mass = np.loadtxt("out_with_bounds/"+prefix+"_0001_aero_size_mass.txt")
    for row in range(0,size_mass.shape[0]):
        for col in range(1,size_mass.shape[1]):
            array_max_size_mass[row,col] = max(array_max_size_mass[row,col],size_mass[row,col])
            array_min_size_mass[row,col] = min(array_min_size_mass[row,col],size_mass[row,col])

# plot time series of number and mass concentrations
partmc_time = np.loadtxt("out_with_bounds/0925_0001_aero_time.txt")
barrel_time = np.loadtxt("ref_aero_time.txt")
(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.plot(barrel_time[:,0], barrel_time[:,1], color='r')
axes.plot(partmc_time[:,0], partmc_time[:,1], color='k')
axes.plot(time, max_time_num, linestyle='--', color='k')
axes.plot(time, min_time_num, linestyle='--', color='k')
axes.set_title("prefactor = %.3f, exponent = %.3f" % (prefactor, exponent))
axes.set_xlabel("Time (s)")
axes.set_ylabel(r"Number concentration ($\mathrm{m}^{-3}$)")
axes.grid()
axes.legend(('Barrel', 'PartMC'))
filename_out = "plot_aero_time_num_with_bounds.pdf"
figure.savefig(filename_out)
(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.plot(barrel_time[:,0], barrel_time[:,2], color='r')
axes.plot(partmc_time[:,0], partmc_time[:,2], color='k')
axes.plot(time, max_time_mass, linestyle='--', color='k')
axes.plot(time, min_time_mass, linestyle='--', color='k')
axes.set_title("prefactor = %.3f, exponent = %.3f" % (prefactor, exponent))
axes.set_xlabel("Time (s)")
axes.set_ylabel(r"Mass concentration (kg $\mathrm{m}^{-3}$)")
axes.grid()
axes.legend(('Barrel', 'PartMC'))
filename_out = "plot_aero_time_mass_with_bounds.pdf"
figure.savefig(filename_out)

# plot number and mass distributions
partmc_num = np.loadtxt("out_with_bounds/0925_0001_aero_size_num.txt")
barrel_num = np.loadtxt("ref_aero_size_num_regrid.txt")
partmc_mass = np.loadtxt("out_with_bounds/0925_0001_aero_size_mass.txt")
barrel_mass = np.loadtxt("ref_aero_size_mass_regrid.txt")
t = 0

if not os.path.exists("plots_with_bounds_num_size"):
   os.mkdir("plots_with_bounds_num_size")
for col in range(1,partmc_num.shape[1]):
     (figure, axes) = mpl_helper.make_fig(colorbar=False)
     axes.semilogx(barrel_num[:,0], barrel_num[:,col], color='r')
     axes.semilogx(partmc_num[:,0], partmc_num[:,col]*math.log(10), color='k')
     axes.semilogx(array_max_size_num[:,0], array_max_size_num[:,col]*math.log(10), linestyle='--', color='k')
     axes.semilogx(array_min_size_num[:,0], array_min_size_num[:,col]*math.log(10), linestyle='--', color='k')
     t = t + 1
     axes.set_title(r"k\_{D} = %.3f, a = %.3f at time %02d" \
     % (prefactor, exponent, t))
     axes.set_xlabel("Dry diameter (m)")
     axes.set_ylabel(r"Number concentration ($\mathrm{m}^{-3}$)")
     axes.grid()
     axes.legend(('Barrel', 'PartMC'))
     filename_out = "plots_with_bounds_num_size/plot_aero_size_num_with_bounds_time%02d.png" %(t)
     #figure.set_dpi(600)
     figure.savefig(filename_out)
#for time in [0,10,20,30]:
#    (figure, axes) = mpl_helper.make_fig(colorbar=False)
#    axes.semilogx(barrel_mass[:,0], barrel_mass[:,time+1])
#    axes.semilogx(partmc_mass[:,0], partmc_mass[:,time+1]*math.log(10), marker='o', mfc='None')
#    axes.set_title("k\_{D} = %.3f, a = %.3f at time %02d" \
#    % (prefactor, exponent, time))
#    axes.set_xlabel("Dry diameter (m)")
#    axes.set_ylabel("Mass concentration (kg $\mathrm{m}^{-3}$)")
#    axes.grid()
#    axes.legend(('Barrel', 'PartMC'), loc='upper left')
#    filename_out = "plot_aero_size_mass_time%02d.pdf" %(time)
#    figure.savefig(filename_out)
