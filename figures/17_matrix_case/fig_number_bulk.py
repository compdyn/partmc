#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

netcdf_dir = "../../local_scenarios/matrix_case/out/"
netcdf_pattern = "brownian_part_0001_(.*).nc"
time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)

matrix_number_a = np.loadtxt('data/set03nmba.txt') 
matrix_number_b = np.loadtxt('data/set03nmbb.txt') 
matrix_so4_a = np.loadtxt('data/set03so4a.txt')
matrix_so4_b = np.loadtxt('data/set03so4b.txt')
matrix_bc = np.loadtxt('data/set03bcar.txt')
matrix_spcm = np.loadtxt('data/set03spcm.txt')

matrix_time = matrix_number_a[:,0]
print matrix_time

matrix_number_akk = matrix_number_a[:,1]
matrix_number_acc = matrix_number_a[:,2]
matrix_number_bc1 = matrix_number_b[:,2]
matrix_number_bc2 = matrix_number_b[:,3]
matrix_number_bc3 = matrix_number_b[:,4]
matrix_number_bcs = matrix_number_b[:,7]
matrix_number_total = matrix_number_akk + matrix_number_acc + matrix_number_bc1 + matrix_number_bc2 + matrix_number_bc3 + matrix_number_bcs

matrix_mass_so4_total = matrix_spcm[:,1]
matrix_mass_bc_total = matrix_spcm[:,2]
matrix_mass_total =  matrix_mass_so4_total + matrix_mass_bc_total
matrix_mass_so4_acc = matrix_so4_a[:,2]
matrix_mass_bc_bc1 = matrix_bc[:,1]
matrix_mass_so4_bcs = matrix_so4_b[:,7]
matrix_mass_bc_bcs = matrix_bc[:,6]

print "number acc t=0 ", matrix_number_acc[0]
print "number bc  t=0 ", matrix_number_bc1[0]
print "mass so4   t=0 ", matrix_mass_so4_total[0]
print "mass bc1   t=0 ", matrix_mass_bc_total[0]

plt.clf()
plt.plot(matrix_time, matrix_number_acc, 'r-', label = 'acc')
plt.plot(matrix_time, matrix_number_bc1, 'k-', label = 'bc1')
plt.plot(matrix_time, matrix_number_bcs, 'b-', label = 'bcs')
plt.plot(matrix_time, matrix_number_total, 'g-', label = 'total')
plt.legend(loc = 'upper right')
plt.xlabel("time ")
plt.ylabel("number concentration in cm^{-3}")
fig = plt.gcf()
fig.savefig("figs/numbers_matrix.pdf")

plt.clf()
plt.plot(matrix_time, matrix_mass_so4_total, 'r-', label = 'so4')
plt.plot(matrix_time, matrix_mass_bc_total, 'k-', label = 'bc')
plt.legend(loc = 'upper right')
plt.xlabel("time ")
plt.ylabel("mass concentration in microgram m^{-3}")
fig = plt.gcf()
fig.savefig("figs/total_mass_matrix.pdf")

array_wc = np.zeros([len(time_filename_list),15])

i_counter = 0
for [time, filename, key] in time_filename_list:
    print time, filename, key
    ncf = scipy.io.netcdf.netcdf_file(filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    env_state = partmc.env_state_t(ncf)
    ncf.close()

    pure_bc = ((particles.masses(include = ["BC"]) > 0) & (particles.masses(include = ["SO4"]) == 0))
    pure_so4 =  ((particles.masses(include = ["SO4"]) > 0) & (particles.masses(include = ["BC"]) == 0))
    with_bc = (particles.masses(include = ["BC"]) > 0)
    with_so4 = (particles.masses(include = ["SO4"]) > 0)
    mixed_bc_so4 = ((particles.masses(include = ["SO4"]) > 0) &  (particles.masses(include = ["BC"]) > 0))
 
    total_number = sum(1/particles.comp_vols)
    pure_so4_number = sum(1/particles.comp_vols[pure_so4])
    pure_bc_number =  sum(1/particles.comp_vols[pure_bc])
    with_so4_number = sum(1/particles.comp_vols[with_so4])
    with_bc_number =  sum(1/particles.comp_vols[with_bc])  
    mixed_bc_so4_number =  sum(1/particles.comp_vols[mixed_bc_so4]) 

    total_mass = sum(particles.masses()/particles.comp_vols)
    bc_total_mass = sum(particles.masses(include = ["BC"])/particles.comp_vols)
    so4_total_mass = sum(particles.masses(include = ["SO4"])/particles.comp_vols)

    bc = particles.masses(include = ["BC"])
    so4 = particles.masses(include = ["SO4"])
    total = particles.masses()
    pure_bc_mass = sum(bc[pure_bc]/particles.comp_vols[pure_bc])
    pure_so4_mass = sum(so4[pure_so4]/particles.comp_vols[pure_so4])
    mixed_bc_so4_mass = sum(total[mixed_bc_so4]/particles.comp_vols[mixed_bc_so4])  
    bc_in_mixed_bc_so4_mass =  sum(bc[mixed_bc_so4]/particles.comp_vols[mixed_bc_so4]) 
    so4_in_mixed_bc_so4_mass =  sum(so4[mixed_bc_so4]/particles.comp_vols[mixed_bc_so4])

    array_wc[i_counter,0]= time / 3600.
    array_wc[i_counter,1]= total_number / 1e6
    array_wc[i_counter,2]= pure_so4_number / 1e6
    array_wc[i_counter,3]= pure_bc_number /1e6
    array_wc[i_counter,4]= with_so4_number /1e6
    array_wc[i_counter,5]= with_bc_number / 1e6
    array_wc[i_counter,6]= mixed_bc_so4_number /1e6
    array_wc[i_counter,7]= total_mass * 1e9
    array_wc[i_counter,8]= bc_total_mass * 1e9
    array_wc[i_counter,9]= so4_total_mass * 1e9
    array_wc[i_counter,10] = pure_bc_mass * 1e9
    array_wc[i_counter,11] = pure_so4_mass * 1e9
    array_wc[i_counter,12] = mixed_bc_so4_mass * 1e9
    array_wc[i_counter,13] = bc_in_mixed_bc_so4_mass * 1e9
    array_wc[i_counter,14] = so4_in_mixed_bc_so4_mass * 1e9
    i_counter += 1

np.savetxt("species_acc_bc1.txt", array_wc)

plt.clf()
plt.plot(array_wc[:,0], array_wc[:,1], 'r-', label = 'total')
plt.plot(array_wc[:,0], array_wc[:,2], 'g-', label = 'pure so4')
plt.plot(array_wc[:,0], array_wc[:,3], 'k-', label = 'pure bc')
#plt.plot(array_wc[:,0], array_wc[:,4], 'm-', label = 'so4-containing')
#plt.plot(array_wc[:,0], array_wc[:,5], 'k-', label = 'bc-containing')
plt.plot(array_wc[:,0], array_wc[:,6], 'b-', label = 'mixed')
plt.plot(matrix_time, matrix_number_total, 'r--', label = 'total')
plt.plot(matrix_time, matrix_number_acc, 'g--', label = 'acc')
plt.plot(matrix_time, matrix_number_bc1, 'k--', label = 'bc1')
plt.plot(matrix_time, matrix_number_bcs, 'b--', label = 'bcs')
plt.axis([0, 24, 0, 2500])
plt.grid(True)
plt.legend(loc = 'upper right')
plt.xlabel("time / h")
plt.ylabel("number concentration in m^{-3}")
fig = plt.gcf()
fig.savefig("figs/numbers.pdf")

plt.clf()
plt.plot(array_wc[:,0], array_wc[:,7], 'b-', label = 'total pmc')
plt.plot(array_wc[:,0], array_wc[:,8], 'k-', label = 'bc total pmc')
plt.plot(array_wc[:,0], array_wc[:,9], 'r-', label = 'so4 total pmc')
plt.plot(matrix_time, matrix_mass_total, 'b--', label = 'total matrix')
plt.plot(matrix_time, matrix_mass_bc_total, 'k--', label = 'bc matrix')
plt.plot(matrix_time, matrix_mass_so4_total, 'r--', label = 'so4 matrix')
plt.axis([0, 24, 0, 6])
plt.grid(True)

plt.legend(loc = 'upper right')
plt.xlabel("time / h ")
plt.ylabel("mass concentration / microgram m^{-3}")
fig = plt.gcf()
fig.savefig("figs/total_masses.pdf")

plt.clf()
plt.plot(array_wc[:,0], array_wc[:,10], 'k-', label = 'bc pure pmc')
plt.plot(array_wc[:,0], array_wc[:,11], 'r-', label = 'so4 pure pmc')
plt.plot(array_wc[:,0], array_wc[:,13], 'g-', label = 'bc mixed pmc')
plt.plot(array_wc[:,0], array_wc[:,14], 'b-', label = 'so4 mixed pmc')
plt.plot(matrix_time, matrix_mass_bc_bc1, 'k--', label = 'bc pure matrix')
plt.plot(matrix_time, matrix_mass_so4_acc, 'r--', label = 'so4 pure matrix')
plt.plot(matrix_time, matrix_mass_bc_bcs, 'g--', label = 'bc mixed matrix')
plt.plot(matrix_time, matrix_mass_so4_bcs, 'b--', label = 'so4 mixed matrix')
plt.axis([0, 24, 0, 6])
plt.grid(True)

plt.legend(loc = 'upper right')
plt.xlabel("time / h ")
plt.ylabel("mass concentration / microgram m^{-3}")
fig = plt.gcf()
fig.savefig("figs/masses_mixing_state.pdf")    

#matrix_mass_so4_total = matrix_spcm[:,1]
#matrix_mass_bc_total = matrix_spcm[:,2]
#matrix_mass_total =  matrix_mass_so4_total + matrix_mass_bc_total
#matrix_mass_so4_acc = matrix_so4_a[:,2]
#matrix_mass_bc_bc1 = matrix_bc[:,1]
#matrix_mass_so4_bcs = matrix_so4_b[:,7]
#matrix_mass_bc_bcs = matrix_bc[:,6]

