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
#for counter in ["10K_flat", "10K_wei-1", "10K_wei-2", "10K_wei-3"]
for counter in ["10K_flat", "10K_wei-1"]
    netcdf_pattern = "urban_plume_wc_%s_0001_(.*).nc"  % counter
    time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)
    time_array = np.zeros([len(time_filename_list)])
    num_avg = np.zeros([len(time_filename_list)])
    mass_avg = np.zeros([len(time_filename_list)])

    array_num = np.zeros([len(time_filename_list),config.i_loop_max])
    array_mass = np.zeros([len(time_filename_list),config.i_loop_max])

    for i_loop in range (0, config.i_loop_max):

        netcdf_pattern = "urban_plume_wc_%s_00%02d_(.*).nc"  % (counter,i_loop+1)
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
            array_mass[i_counter,i_loop]= total_dry_mass * 1e9
            i_counter += 1

    num_avg = np.average(array_num, axis = 1)
    num_std = np.std(array_num, axis = 1)
    mass_avg = np.average(array_mass, axis = 1)
    mass_std = np.std(array_mass, axis = 1)

    print 'average of number std ', np.average(num_std)
    print 'average of mass std ', np.average(mass_std)

    f1 = "data/time_array_%s.txt" % counter
    f2 = "data/array_number_%s.txt" % counter
    f3 = "data/array_mass_%s.txt" % counter
    f4 = "data/average_number_%s.txt" % counter
    f5 = "data/average_mass_%s.txt" % counter
    f6 = "data/std_number_%s.txt" % counter
    f7 = "data/std_mass_%s.txt" % counter

    np.savetxt(f1, time_array)
    np.savetxt(f2, array_num)
    np.savetxt(f3, array_mass)
    np.savetxt(f4, num_avg)
    np.savetxt(f5, mass_avg)
    np.savetxt(f6, num_std)
    np.savetxt(f7, mass_std)




    



