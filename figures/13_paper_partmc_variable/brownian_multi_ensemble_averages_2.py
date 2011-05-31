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

i_ens_max = config.i_ens_max
i_loop_max = config.i_loop_max
i_bin = 100

sect_array_num_all = np.loadtxt("/home/nriemer/subversion/partmc/branches/nriemer/local_scenarios/brownian_test_paper2/brownian_sect_size_num.txt")
sect_array_mass_all = np.loadtxt("/home/nriemer/subversion/partmc/branches/nriemer/local_scenarios/brownian_test_paper2/brownian_sect_size_mass.txt")

sect_array_num = np.log(10) * sect_array_num_all[:,12]
sect_array_mass = np.log(10) * sect_array_mass_all[:,12]
print "sect_array_num ", sect_array_num.shape

sect_num = sect_array_num.reshape((100,10)).mean(axis=1)
sect_mass = sect_array_mass.reshape((100,10)).mean(axis=1)

norm_num = np.zeros([i_ens_max,i_loop_max])
norm_mass = np.zeros([i_ens_max,i_loop_max])

norm_num_avg = np.zeros([i_loop_max])
norm_mass_avg = np.zeros([i_loop_max])

num_avg = np.zeros([i_ens_max,i_loop_max,i_bin])
num_std = np.zeros([i_ens_max,i_loop_max,i_bin])
mass_avg = np.zeros([i_ens_max,i_loop_max,i_bin])
mass_std = np.zeros([i_ens_max,i_loop_max,i_bin])

f1 = "data/ensemble_size_dist_num_brownian_600s.txt"
f2 = "data/ensemble_size_dist_mass_brownian_600s.txt" 

array_num = np.loadtxt(f1)
array_mass = np.loadtxt(f2)

for i_ens in range(0, i_ens_max):
    for i_loop in range(0, i_loop_max):
        print i_ens, i_loop
        array_index_start = i_ens * 100
        array_index_end = i_loop + 1 + i_ens * 100
        print "array_index start and end ", array_index_start, array_index_end
        num_avg[i_ens,i_loop,:] = np.average(array_num[array_index_start:array_index_end,:], axis = 0)
        num_std[i_ens,i_loop,:] = np.std(array_num[array_index_start:array_index_end,:], axis = 0)
        
        mass_avg[i_ens,i_loop,:] = np.average(array_mass[array_index_start:array_index_end,:], axis = 0)
        mass_std[i_ens,i_loop,:] = np.std(array_mass[array_index_start:array_index_end,:], axis = 0)
print "num_avg ", num_avg.shape

for i_ens in range(0,i_ens_max):
    for i_loop in range(0,i_loop_max):
        norm_num[i_ens, i_loop] = np.linalg.norm(num_avg[i_ens, i_loop,:] - sect_num[:])
        norm_mass[i_ens, i_loop] = np.linalg.norm(mass_avg[i_ens, i_loop,:] - sect_mass[:])

norm_num_avg = np.average(norm_num, axis = 0)
norm_num_std = np.std(norm_num, axis = 0)
norm_mass_avg = np.average(norm_mass, axis = 0)
norm_mass_std = np.std(norm_mass, axis = 0)

print "norm_num_avg ", norm_num_avg.shape

f1 = "data/avg_ensemble_norm_num_brownian_600s.txt" 
f2 = "data/avg_ensemble_norm_mass_brownian_600s.txt"
f3 = "data/std_ensemble_norm_num_brownian_600s.txt" 
f4 = "data/std_ensemble_norm_mass_brownian_600s.txt"

np.savetxt(f1, norm_num_avg)
np.savetxt(f2, norm_mass_avg)
np.savetxt(f3, norm_num_std)
np.savetxt(f4, norm_mass_std)


    



