#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc
import config

i_loop_max = config.i_loop_max

netcdf_dir = "../../scenarios/5_weighted/out"
for counter in ["10K_wei\\+1", "10K_flat", "10K_wei-1", "10K_wei-2", "10K_wei-3", "10K_wei-4","1K_wei\\+1", "1K_flat", "1K_wei-1", "1K_wei-2", "1K_wei-3", "1K_wei-4"]:
    netcdf_pattern = "urban_plume_wc_%s_0001_(.*).nc"  % counter
    time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)
    time_array = np.zeros([len(time_filename_list)])

    num_avg = np.zeros([len(time_filename_list), i_loop_max])
    mass_avg = np.zeros([len(time_filename_list), i_loop_max])

    num_std = np.zeros([len(time_filename_list), i_loop_max])
    mass_std = np.zeros([len(time_filename_list), i_loop_max])

    array_num = np.zeros([len(time_filename_list), i_loop_max])
    array_mass = np.zeros([len(time_filename_list), i_loop_max])

    num_avg_overall = np.zeros([i_loop_max])
    mass_avg_overall = np.zeros([i_loop_max])

    num_std_overall = np.zeros([i_loop_max])
    mass_std_overall = np.zeros([i_loop_max])

    for i_loop in range (0, i_loop_max):

        netcdf_pattern = "urban_plume_wc_%s_0%03d_(.*).nc"  % (counter,i_loop+1)
        print netcdf_pattern
        time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)

        i_counter = 0
        for [time, filename, key] in time_filename_list:
            print time, filename, key
            ncf = scipy.io.netcdf.netcdf_file(filename, 'r')
            particles = partmc.aero_particle_array_t(ncf)
            env_state = partmc.env_state_t(ncf)
            ncf.close()

            total_number = sum(1/particles.comp_vols)
            total_dry_mass = sum(particles.masses(exclude = ["H2O"])/particles.comp_vols)
            time_array[i_counter]= time / 3600.
            array_num[i_counter,i_loop]= total_number
            array_mass[i_counter,i_loop]= total_dry_mass * 1e9
            i_counter += 1

    for i_ensemble in range(1, i_loop_max):
        for i_loop in range(0,i_ensemble):
            num_avg[:,i_ensemble] = np.average(array_num[:,:i_loop], axis = 1)
            num_std[:,i_ensemble] = np.std(array_num[:,:i_loop], axis = 1)

            mass_avg[:,i_ensemble] = np.average(array_mass[:,:i_loop], axis = 1)
            mass_std[:,i_ensemble] = np.std(array_mass[:,:i_loop], axis = 1)/mass_avg[:,i_ensemble]

            
    num_avg_overall = np.average(num_avg, axis = 0)
    num_std_overall = np.average(num_std, axis = 0)
    mass_avg_overall = np.average(mass_avg, axis = 0)
    mass_std_overall = np.average(mass_std, axis = 0)

    print 'average of number', num_avg_overall
    print 'average of mass', mass_avg_overall

    f1 = "data/ensemble_number_%s.txt" % counter
    f2 = "data/ensemble_mass_%s.txt" % counter
    f3 = "data/ensemble_number_std_%s.txt" % counter
    f4 = "data/ensemble_mass_std_%s.txt" % counter

    np.savetxt(f1, num_avg_overall)
    np.savetxt(f2, mass_avg_overall)
    np.savetxt(f3, num_std_overall)
    np.savetxt(f4, mass_std_overall)




    



