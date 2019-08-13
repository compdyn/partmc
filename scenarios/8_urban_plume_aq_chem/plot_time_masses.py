#!/usr/bin/env python

import scipy.io
import sys
sys.path.append("../../tool")
import numpy as np
import mpl_helper
import matplotlib
import matplotlib.cm
import matplotlib.ticker
import matplotlib.pyplot as plt
import partmc


bc = np.zeros([601])
oc = np.zeros([601])
so4 = np.zeros([601])
no3 = np.zeros([601])
nh4 = np.zeros([601])
soa = np.zeros([601])
h2o = np.zeros([601])
num_conc = np.zeros([601])
cloud_drop = np.zeros([601])
time = np.zeros([601])
conv_fac = np.zeros([601])


in_dir = "out/"
in_file_pattern = "urban_plume_aq_chem_mono_0001_.*.nc"
env_state_history = partmc.read_history(partmc.env_state_t, in_dir, in_file_pattern)
pressure = np.array([env_state_history[i][1].pressure for i in range(len(env_state_history))])
temp = np.array([env_state_history[i][1].temperature for i in range(len(env_state_history))])

print(pressure)
print(temp)

R_d = 287.058
conv_fac = pressure / R_d / temp
print(conv_fac)

for counter in range(0, 601):
    print("counter = ",  counter)
    time[counter] = counter
    print("time = ", time[counter])
    
    in_filename = "out/urban_plume_aq_chem_mono_0001_0000%04d.nc" % (counter+1)
    print("in file ", in_filename)
    ncf = scipy.io.netcdf.netcdf_file(in_filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    ncf.close()

    wet_diameters = particles.diameters() * 1e6
    num_conc[counter] = sum(particles.num_concs) / conv_fac[counter]

    activated = (wet_diameters > 2e0)

    cloud_drop[counter] = sum(particles.num_concs[activated]) / conv_fac[counter]
    
    bc[counter] = sum(particles.masses(include = ["BC"])*particles.num_concs) / conv_fac[counter]
    h2o[counter] = sum(particles.masses(include = ["H2O"])*particles.num_concs) / conv_fac[counter]
    oc[counter] = sum(particles.masses(include = ["OC"])*particles.num_concs) / conv_fac[counter]
    so4[counter] = sum(particles.masses(include = ["SO4"])*particles.num_concs) / conv_fac[counter]
    nh4[counter] = sum(particles.masses(include = ["NH4"])*particles.num_concs) / conv_fac[counter]
    no3[counter] = sum(particles.masses(include = ["NO3"])*particles.num_concs) / conv_fac[counter]
    soa[counter] = sum(particles.masses(include = ["API1", "API2", "LIM1", "LIM2", "ARO1", "ARO2", "OLE1", "ALK1"])*particles.num_concs) /conv_fac[counter]

    
    
(figure, axes) = mpl_helper.make_fig(figure_width=5,
                                     top_margin=0.5, bottom_margin=0.45,
                                     left_margin=0.65, right_margin=1.2)

axes.set_xscale("linear")
axes.set_yscale("linear")
axes.plot(time/60., h2o*1e3, label = 'H2O')

axes.set_ylabel(r"water mass conc. in $\rm g\, kg^{-3}$")
axes.set_xlabel(r"time in minutes")
#axes.set_xlim(0,50)
#axes.set_xticks([0, 10, 20, 30, 40, 50])
axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
axes.grid(True)
figure.savefig("figs/time_h2o_mono.pdf")
plt.close()

(figure, axes) = mpl_helper.make_fig(figure_width=5,
                                     top_margin=0.5, bottom_margin=0.45,
                                     left_margin=0.65, right_margin=1.2)

axes.set_xscale("linear")
axes.set_yscale("linear")
axes.plot(time/60, num_conc/1e6, label = 'num. conc')
axes.plot(time/60, cloud_drop/1e6, label = 'cloud drop. conc.')
axes.set_ylabel(r"number conc. in $\rm kg^{-3}$")
axes.set_xlabel(r"time in minutes")
#axes.set_xlim(0,50)
#axes.set_xticks([0, 10, 20, 30, 40, 50])
axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
axes.grid(True)
figure.savefig("figs/time_number_conc_mono.pdf")
plt.close()


(figure, axes) = mpl_helper.make_fig(figure_width=5,
                                     top_margin=0.5, bottom_margin=0.45,
                                     left_margin=0.65, right_margin=1.2)

axes.set_xscale("linear")
axes.set_yscale("linear")
axes.plot(time/60, bc*1e9, label = 'BC')
axes.plot(time/60, oc*1e9, label = 'OC')
axes.plot(time/60, so4*1e9, label = 'SO4')
axes.plot(time/60, no3*1e9, label = 'NO3')
axes.plot(time/60, soa*1e9, label = 'SOA')
axes.plot(time/60, nh4*1e9, label = 'NH4')

axes.set_ylabel(r"mass conc. in $\rm \mu g\, kg^{-3}$")
axes.set_xlabel(r"time in minutes")
#axes.set_xlim(0,50)
#axes.set_xticks([0, 10, 20, 30, 40, 50])
axes.legend(loc='center left', bbox_to_anchor=(1, 0.5))
axes.grid(True)
figure.savefig("figs/time_masses_mono.pdf")
plt.close()

   


