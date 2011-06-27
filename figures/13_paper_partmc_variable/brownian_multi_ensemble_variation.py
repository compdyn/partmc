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
i_ens_max = config.i_ens_max

netcdf_dir = "/home/nriemer/subversion/partmc/branches/nriemer/local_scenarios/brownian_test_paper2/out/"

array_num_init = np.zeros([i_loop_max*i_ens_max])
array_mass_init = np.zeros([i_loop_max*i_ens_max])

array_num_end = np.zeros([i_loop_max*i_ens_max])
array_mass_end = np.zeros([i_loop_max*i_ens_max])

for i_loop in range (0, i_ens_max*i_loop_max):
    filename = "brownian_part_%05d_00000001.nc"  % (i_loop+1)
    print filename
    ncf = scipy.io.netcdf.netcdf_file(netcdf_dir+filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    env_state = partmc.env_state_t(ncf)
    ncf.close()

    total_number = sum(1/particles.comp_vols)
    total_dry_mass = sum(particles.masses()/particles.comp_vols)

    array_num_init[i_loop] = total_number
    array_mass_init[i_loop]= total_dry_mass
 
    filename = "brownian_part_%05d_00000002.nc"  % (i_loop+1)
    print filename
    ncf = scipy.io.netcdf.netcdf_file(netcdf_dir+filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    env_state = partmc.env_state_t(ncf)
    ncf.close()

    total_number = sum(1/particles.comp_vols)
    total_dry_mass = sum(particles.masses()/particles.comp_vols)

    array_num_end[i_loop] = total_number
    array_mass_end[i_loop]= total_dry_mass

    
delta_num = array_num_init - array_num_end
delta_mass = array_mass_init - array_mass_end

num_avg_init = np.average(array_num_init, axis = 0)
num_std_init = np.std(array_num_init, axis = 0)
mass_avg_init = np.average(array_mass_init, axis = 0)
mass_std_init = np.std(array_mass_init, axis = 0)

delta_num_avg = np.average(delta_num, axis = 0)
delta_num_std = np.std(delta_num, axis = 0)
delta_mass_avg = np.average(delta_mass, axis = 0)
delta_mass_std = np.std(delta_mass, axis = 0)

print
print "number avg init ", num_avg_init
print "number std init ", num_std_init
print "delta_num avg ", delta_num_avg
print "delta_num std ", delta_num_std
print
print "mass avg init ", mass_avg_init
print "mass std init ", mass_std_init
print "delta_mass avg ", delta_mass_avg
print "delta_mass std ", delta_mass_std
print
print "ratio std number ", delta_num_std / num_std_init
print "ratio std mass ", delta_mass_std / mass_std_init



    



