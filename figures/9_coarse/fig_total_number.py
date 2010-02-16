#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc
import config

netcdf_dir = "../../scenarios/5_weighted/out"
netcdf_pattern = "urban_plume_wc_10K_flat_0001_(.*).nc"
time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)
time_array = np.zeros([len(time_filename_list)])
num_avg = np.zeros([len(time_filename_list)])
mass_avg = np.zeros([len(time_filename_list)])


for i_loop in range (0, config.i_loop_max):

    netcdf_pattern = "urban_plume_wc_10K_flat_00%02d_(.*).nc"  % (i_loop+1)
    print netcdf_pattern
    time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)
    
    array_num = np.zeros([len(time_filename_list),config.i_loop_max])
    array_mass = np.zeros([len(time_filename_list),config.i_loop_max])

    i_counter = 0
    for [time, filename, key] in time_filename_list:
        print time, filename, key
        ncf = Scientific.IO.NetCDF.NetCDFFile(filename)
        particles = partmc.aero_particle_array_t(ncf)
        env_state = partmc.env_state_t(ncf)
        ncf.close()

        total_number = sum(1/particles.comp_vols)
        total_dry_mass = sum(particles.masses(exclude = ["H2O"])/particles.comp_vols)
        time_array[i_counter]= time / 3600.
        array_num[i_counter,i_loop]= total_number
        print 'total_number ', i_counter, i_loop, total_number, array_num[i_counter,i_loop]
        array_mass[i_counter,i_loop]= total_dry_mass * 1e9
        i_counter += 1

print array_num.shape, array_mass.shape, time_array.shape
num_avg = np.average(array_num, axis = 1)
mass_avg = np.average(array_mass, axis = 1)

print 'numbers ', array_num[:,0]
plt.clf()
for i_loop in range(0,config.i_loop_max):
    plt.plot(time_array[:], array_num[:,i_loop], 'k')
plt.plot(time_array[:], num_avg[:], 'r')
plt.xlabel("time ")
plt.ylabel("number concentration in m^{-3}")
plt.title("10K flat")
fig = plt.gcf()
fig.savefig("figs/number_10K_flat.pdf")

plt.clf()
for i_loop in range(0,config.i_loop_max):
    plt.plot(time_array[:], array_mass[:,i_loop], 'k')
plt.plot(time_array[:], mass_avg[:], 'r')
plt.xlabel("time ")
plt.ylabel("mass concentration in \mu g m^{-3}")
plt.title("10K flat")
fig = plt.gcf()
fig.savefig("figs/mass_10K_flat.pdf")




    



