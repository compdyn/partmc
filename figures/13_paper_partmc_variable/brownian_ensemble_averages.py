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

sect_array_num_all = np.loadtxt("/Users/nriemer/subversion/partmc/trunk/local_scenarios/brownian_test_paper/out_60s/brownian_sect_size_num.txt")
sect_array_mass_all = np.loadtxt("/Users/nriemer/subversion/partmc/trunk/local_scenarios/brownian_test_paper/out_60s/brownian_sect_size_mass.txt")

sect_array_num = np.log(10) * sect_array_num_all[:,25]
sect_array_mass = np.log(10) * sect_array_mass_all[:,25]
print "sect_array_num ", sect_array_num.shape

sect_num = sect_array_num.reshape((100,10)).mean(axis=1)
sect_mass = sect_array_mass.reshape((100,10)).mean(axis=1)

print "sect_num ", sect_num.shape

norm_num = np.zeros([i_loop_max])
norm_mass = np.zeros([i_loop_max])

x_axis = partmc.log_grid(min=1e-10,max=1e-4,n_bin=100)
x_centers = x_axis.centers()

netcdf_pattern = "brownian_part_0001_(.*).nc"

time_filename_list = partmc.get_time_filename_list(config.netcdf_dir, netcdf_pattern)
time_array = np.zeros([len(time_filename_list)])

num_avg = np.zeros([len(time_filename_list), i_loop_max, len(x_centers)])
mass_avg = np.zeros([len(time_filename_list), i_loop_max, len(x_centers)])

num_std = np.zeros([len(time_filename_list), i_loop_max, len(x_centers)])
mass_std = np.zeros([len(time_filename_list), i_loop_max, len(x_centers)])

array_num = np.zeros([len(time_filename_list), i_loop_max, len(x_centers)])
array_mass = np.zeros([len(time_filename_list), i_loop_max, len(x_centers)])

num_avg_overall = np.zeros([i_loop_max, len(x_centers)])
mass_avg_overall = np.zeros([i_loop_max, len(x_centers)])

num_std_overall = np.zeros([i_loop_max, len(x_centers)])
mass_std_overall = np.zeros([i_loop_max, len(x_centers)])

for i_loop in range (0, i_loop_max):

    netcdf_pattern = "brownian_part_0%03d_(.*).nc"  % (i_loop+1)
    print netcdf_pattern
    time_filename_list = partmc.get_time_filename_list(config.netcdf_dir, netcdf_pattern)

    i_counter = 0
    for [time, filename, key] in time_filename_list:
        print time, filename, key
        ncf = scipy.io.netcdf.netcdf_file(filename, 'r')
        particles = partmc.aero_particle_array_t(ncf)
        env_state = partmc.env_state_t(ncf)
        ncf.close()

        wet_diameters = particles.diameters()
        hist = partmc.histogram_1d(wet_diameters, x_axis, weights = 1 / particles.comp_vols)

#        total_number = sum(1/particles.comp_vols)
#        total_mass = sum(particles.masses()/particles.comp_vols)
        time_array[i_counter]= time / 3600.
        array_num[i_counter,i_loop,:]= hist

        hist = partmc.histogram_1d(wet_diameters, x_axis, weights = particles.masses() / particles.comp_vols)
        array_mass[i_counter,i_loop,:]= hist
        i_counter += 1

print "array_num ", array_num.shape

for i_ensemble in range(0, i_loop_max):
    print "i_ensemble ", i_ensemble
    num_avg[:,i_ensemble,:] = np.average(array_num[:,:i_ensemble+1,:], axis = 1)
    num_std[:,i_ensemble,:] = np.std(array_num[:,:i_ensemble+1,:], axis = 1)

    mass_avg[:,i_ensemble,:] = np.average(array_mass[:,:i_ensemble+1,:], axis = 1)
    mass_std[:,i_ensemble,:] = np.std(array_mass[:,:i_ensemble+1,:], axis = 1)

print "num_avg ", num_avg.shape

for i_ensemble in range(0,i_loop_max):
    norm_num[i_ensemble] = np.linalg.norm(num_avg[24,i_ensemble,:] - sect_num[:])
    norm_mass[i_ensemble] = np.linalg.norm(mass_avg[24,i_ensemble,:] - sect_mass[:])

print "norms ", norm_num, norm_mass
f1 = "data/ensemble_norm_num_brownian_1200s.txt"
f2 = "data/ensemble_norm_mass_brownian_1200s.txt"

np.savetxt(f1, norm_num)
np.savetxt(f2, norm_mass)



    



