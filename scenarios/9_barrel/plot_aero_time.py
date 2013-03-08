#!/usr/bin/env python
import sys
sys.path.append('../../tool/')
import subprocess
import os
import math
import numpy as np
import mpl_helper
import matplotlib.pyplot as plt

shape_opt = raw_input("Enter s for spherical, f for fractal:")

partmc_time = np.loadtxt("out_0212/barrel_wc_nummass_source_0001_aero_time.txt")
barrel_time = np.loadtxt("ref_0212/ref_aero_time.txt")
(figure, axes) = mpl_helper.make_fig(colorbar=False)
axes.plot(barrel_time[:,0], barrel_time[:,1], color='r')
axes.plot(partmc_time[:,0], partmc_time[:,1], color='k')
axes.set_xlabel("Time (s)")
axes.set_ylabel(r"Number concentration ($\mathrm{m}^{-3}$)")
axes.grid()
axes.legend(('Barrel', 'PartMC'))
filename_out = "aero_num_time.pdf"
figure.savefig(filename_out)
if shape_opt == 's':
   (figure, axes) = mpl_helper.make_fig(colorbar=False)
   axes.plot(barrel_time[:,0], barrel_time[:,2], color='r')
   axes.plot(partmc_time[:,0], partmc_time[:,2], color='k')
   axes.set_xlabel("Time (s)")
   axes.set_ylabel(r"Mass concentration (kg $\mathrm{m}^{-3}$)")
   axes.grid()
   axes.legend(('Barrel', 'PartMC'))
   filename_out = "aero_mass_time.pdf"
   figure.savefig(filename_out)
if shape_opt == 'f':
   rho = 1760 # density in kgm-3
   partmc_num = np.loadtxt("out_0212/barrel_wc_nummass_source_0001_aero_size_num.txt")
   barrel_num = np.loadtxt("ref_0212/ref_aero_size_num_regrid.txt")
   list_partmc_mass = []
   list_barrel_mass = []
   for col in range(1,partmc_num.shape[1]):
       partmc_mass = partmc_num[:,col] * math.pi / 6. * rho * partmc_num[:,0]**3 * math.log(10)
       barrel_mass = barrel_num[:,col] * math.pi / 6. * rho * barrel_num[:,0]**3
       list_partmc_mass.append(sum(partmc_mass * (math.log10(partmc_num[1,0]) - math.log10(partmc_num[0,0]))))
       list_barrel_mass.append(sum(barrel_mass * (math.log10(barrel_num[1,0]) - math.log10(barrel_num[0,0]))))
   (figure, axes) = mpl_helper.make_fig(colorbar=False)
   axes.plot(barrel_time[:,0], list_barrel_mass, color='r')
   axes.plot(partmc_time[:,0], list_partmc_mass, color='k')
   axes.set_xlabel("Time (s)")
   axes.ticklabel_format(style='sci', scilimits=(0,0), axis='y')
   axes.set_ylabel(r"Mass concentration (kg $\mathrm{m}^{-3}$)")
   axes.grid()
   axes.legend(('Barrel', 'PartMC'))
   filename_out = "aero_mass_time.pdf"
   figure.savefig(filename_out)
