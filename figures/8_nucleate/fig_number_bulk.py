#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc
const = pmc_data_nc.load_constants("../../src/constants.f90")

netcdf_dir = "../../urban_plume_nucleate/out/"
netcdf_pattern = "urban_plume_wc_0001_(.*).nc"
time_filename_list = pmc_data_nc.get_time_filename_list(netcdf_dir, netcdf_pattern)

array_wc = np.zeros([len(time_filename_list),9])

i_counter = 0
for [time, filename, key] in time_filename_list:
    print time, filename, key
    ncf = Scientific.IO.NetCDF.NetCDFFile(filename)
    particles = pmc_data_nc.aero_particle_array_t(ncf)
    env_state = pmc_data_nc.env_state_t(ncf)
    ncf.close()

    total_number = sum(1/particles.comp_vol)
    total_dry_mass = sum(particles.mass(exclude = ["H2O"])/particles.comp_vol)
    total_mass = sum(particles.mass()/particles.comp_vol)
    bc = sum(particles.mass(include = ["BC"])/particles.comp_vol)
    oc = sum(particles.mass(include = ["OC"])/particles.comp_vol)
    so4 = sum(particles.mass(include = ["SO4"])/particles.comp_vol)
    nh4 = sum(particles.mass(include = ["NH4"])/particles.comp_vol)
    no3 = sum(particles.mass(include = ["NO3"])/particles.comp_vol)

    array_wc[i_counter,0]= time / 3600.
    array_wc[i_counter,1]= total_number
    array_wc[i_counter,2]= total_dry_mass
    array_wc[i_counter,3]= total_mass
    array_wc[i_counter,4]= bc
    array_wc[i_counter,5]= oc
    array_wc[i_counter,6]= so4
    array_wc[i_counter,7]= nh4
    array_wc[i_counter,8]= no3
    
    i_counter += 1

plt.clf()
plt.plot(array_wc[:,0], array_wc[:,1], 'r-', label = 'wc')
plt.legend(loc = 'upper right')
plt.xlabel("time ")
plt.ylabel("number concentration in m^{-3}")
fig = plt.gcf()
fig.savefig("figs/number_wc.pdf")

plt.clf()
plt.plot(array_wc[:,0], array_wc[:,2], 'r-', label = 'dry wc')
plt.plot(array_wc[:,0], array_wc[:,3], 'b-', label = 'total wc')
plt.legend(loc = 'upper right')
plt.xlabel("time ")
plt.ylabel("mass concentration in kg m^{-3}")
fig = plt.gcf()
fig.savefig("figs/mass_wc.pdf")

plt.clf()
plt.plot(array_wc[:,0], array_wc[:,4], 'r-', label = 'bc wc')
plt.plot(array_wc[:,0], array_wc[:,5], 'b-', label = 'oc wc')
plt.plot(array_wc[:,0], array_wc[:,6], 'g-', label = 'so4 wc')
plt.plot(array_wc[:,0], array_wc[:,7], 'k-', label = 'nh4 wc')
plt.plot(array_wc[:,0], array_wc[:,8], 'm-', label = 'no3 wc')
plt.legend(loc = 'upper right')
plt.xlabel("time ")
plt.ylabel("mass concentration in kg m^{-3}")
fig = plt.gcf()
fig.savefig("figs/mass_species_wc.pdf")




    



