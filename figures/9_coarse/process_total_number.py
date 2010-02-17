#!/usr/bin/env python

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

array_num = np.zeros([len(time_filename_list),config.i_loop_max])
array_mass = np.zeros([len(time_filename_list),config.i_loop_max])

for i_loop in range (0, config.i_loop_max):

    netcdf_pattern = "urban_plume_wc_10K_wei-3_00%02d_(.*).nc"  % (i_loop+1)
    print netcdf_pattern
    time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)
    
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
num_std = np.std(array_num, axis = 1)
mass_avg = np.average(array_mass, axis = 1)
mass_std = np.std(array_mass, axis = 1)

print 'average of number std ', np.average(num_std)
print 'average of mass std ', np.average(mass_std)

np.savetxt("data/time_array.txt", time_array)
np.savetxt("data/array_number.txt", array_num)
np.savetxt("data/array_mass.txt", array_mass)
np.savetxt("data/average_number.txt", num_avg)
np.savetxt("data/average_mass.txt", mass_avg)
np.savetxt("data/std_number.txt", num_std)
np.savetxt("data/std_mass.txt", mass_std)




    



