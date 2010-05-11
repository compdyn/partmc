#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

netcdf_dir = "../../scenarios/5_weighted/out"
netcdf_pattern = "urban_plume_wc_10K_flat_0001_(.*).nc"
time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)

array_wc = np.zeros([len(time_filename_list),10])
array_nc = np.zeros([len(time_filename_list),10])

i_counter = 0
for [time, filename, key] in time_filename_list:
    print time, filename, key
    ncf = scipy.io.netcdf.netcdf_file(filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    env_state = partmc.env_state_t(ncf)
    ncf.close()

    total_number = sum(1/particles.comp_vols)
    total_dry_mass = sum(particles.masses(exclude = ["H2O"])/particles.comp_vols)
    total_mass = sum(particles.masses()/particles.comp_vols)
    bc = sum(particles.masses(include = ["BC"])/particles.comp_vols)
    oc = sum(particles.masses(include = ["OC"])/particles.comp_vols)
    so4 = sum(particles.masses(include = ["SO4"])/particles.comp_vols)
    nh4 = sum(particles.masses(include = ["NH4"])/particles.comp_vols)
    no3 = sum(particles.masses(include = ["NO3"])/particles.comp_vols)
    oin = sum(particles.masses(include = ["OIN"])/particles.comp_vols)

    array_wc[i_counter,0]= time / 3600.
    array_wc[i_counter,1]= total_number
    array_wc[i_counter,2]= total_dry_mass
    array_wc[i_counter,3]= total_mass
    array_wc[i_counter,4]= bc
    array_wc[i_counter,5]= oc
    array_wc[i_counter,6]= so4
    array_wc[i_counter,7]= nh4
    array_wc[i_counter,8]= no3
    array_wc[i_counter,9]= oin
    
    i_counter += 1

#netcdf_dir = "../../scenarios/5_weighted/out/"
#netcdf_pattern = "urban_plume_nc_0001_(.*).nc"
#time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)

#i_counter = 0
#for [time, filename, key] in time_filename_list:
#    print time, filename, key
#    ncf = scipy.io.netcdf.netcdf_file(filename, 'r')
#    particles = partmc.aero_particle_array_t(ncf)
#    env_state = partmc.env_state_t(ncf)
#    ncf.close()
#
#    total_number = sum(1/particles.comp_vols)
#    total_dry_mass = sum(particles.masses(exclude = ["H2O"])/particles.comp_vols)
#    total_mass = sum(particles.masses()/particles.comp_vols)
#    bc = sum(particles.masses(include = ["BC"])/particles.comp_vols)
#    oc = sum(particles.masses(include = ["OC"])/particles.comp_vols)
#    so4 = sum(particles.masses(include = ["SO4"])/particles.comp_vols)
#    nh4 = sum(particles.masses(include = ["NH4"])/particles.comp_vols)
#    no3 = sum(particles.masses(include = ["NO3"])/particles.comp_vols)
#    oin = sum(particles.masses(include = ["OIN"])/particles.comp_vols)
#
#    array_nc[i_counter,0]= time / 3600.
#    array_nc[i_counter,1]= total_number
#    array_nc[i_counter,2]= total_dry_mass
#    array_nc[i_counter,3]= total_mass
#    array_nc[i_counter,4]= bc
#    array_nc[i_counter,5]= oc
#    array_nc[i_counter,6]= so4
#    array_nc[i_counter,7]= nh4
#    array_nc[i_counter,8]= no3
#    array_nc[i_counter,9]= oin
#    i_counter += 1

plt.clf()
plt.plot(array_wc[:,0], array_wc[:,1], 'r-', label = 'wc')
#plt.plot(array_nc[:,0], array_nc[:,1], 'r--', label = 'nc')
plt.legend(loc = 'upper right')
plt.xlabel("time ")
plt.ylabel("number concentration in m^{-3}")
fig = plt.gcf()
fig.savefig("figs/number_coarse_w.pdf")

plt.clf()
plt.plot(array_wc[:,0], array_wc[:,2], 'r-', label = 'dry wc')
plt.plot(array_wc[:,0], array_wc[:,3], 'b-', label = 'total wc')
#plt.plot(array_nc[:,0], array_nc[:,2], 'r--', label = 'dry nc')
#plt.plot(array_nc[:,0], array_nc[:,3], 'b--', label = 'total nc')
plt.legend(loc = 'upper right')
plt.xlabel("time ")
plt.ylabel("mass concentration in kg m^{-3}")
fig = plt.gcf()
fig.savefig("figs/mass_coarse_w.pdf")

plt.clf()
plt.plot(array_wc[:,0], array_wc[:,4], 'r-', label = 'bc wc')
plt.plot(array_wc[:,0], array_wc[:,5], 'b-', label = 'oc wc')
plt.plot(array_wc[:,0], array_wc[:,6], 'g-', label = 'so4 wc')
plt.plot(array_wc[:,0], array_wc[:,7], 'k-', label = 'nh4 wc')
plt.plot(array_wc[:,0], array_wc[:,8], 'm-', label = 'no3 wc')
plt.plot(array_wc[:,0], array_wc[:,9], 'm-', label = 'oin wc')
plt.legend(loc = 'upper right')
plt.xlabel("time ")
plt.ylabel("mass concentration in kg m^{-3}")
fig = plt.gcf()
fig.savefig("figs/mass_species_wc_coarse_w.pdf")

#plt.clf()
#plt.plot(array_nc[:,0], array_nc[:,4], 'r-', label = 'bc nc')
#plt.plot(array_nc[:,0], array_nc[:,5], 'b-', label = 'oc nc')
#plt.plot(array_nc[:,0], array_nc[:,6], 'g-', label = 'so4 nc')
#plt.plot(array_nc[:,0], array_nc[:,7], 'k-', label = 'nh4 nc')
#plt.plot(array_nc[:,0], array_nc[:,8], 'm-', label = 'no3 nc')
#plt.plot(array_nc[:,0], array_nc[:,9], 'm-', label = 'oin nc')
#plt.legend(loc = 'upper right')
#plt.xlabel("time ")
#plt.ylabel("mass concentration in kg m^{-3}")
#fig = plt.gcf()
#fig.savefig("figs/mass_species_nc_coarse_w.pdf")



    



