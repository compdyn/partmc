#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc
const = partmc.constants_t("../../src/constants.f90")

netcdf_dir = "../../scenarios/5_weighted/out"
netcdf_pattern1 = "urban_plume_wc_0001_(.*).nc"
netcdf_pattern2 = "urban_plume_wc_0002_(.*).nc"
netcdf_pattern3 = "urban_plume_wc_0003_(.*).nc"

time_filename_list1 = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern1)
time_filename_list2 = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern2)
time_filename_list3 = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern3)

time_array = np.zeros([len(time_filename_list1)])

array_num = np.zeros([len(time_filename_list1),3])
array_mass = np.zeros([len(time_filename_list1),3])

num_avg = np.zeros([len(time_filename_list1)])
mass_avg = np.zeros([len(time_filename_list1)])

i_counter = 0
for [time, filename, key] in time_filename_list1:
    print time, filename, key
    ncf = Scientific.IO.NetCDF.NetCDFFile(filename)
    particles = partmc.aero_particle_array_t(ncf)
    env_state = partmc.env_state_t(ncf)
    ncf.close()

    total_number = sum(1/particles.comp_vols)
    total_dry_mass = sum(particles.mass(exclude = ["H2O"])/particles.comp_vols)
    time_array[i_counter]= time / 3600.
    array_num[i_counter,0]= total_number
    array_mass[i_counter,0]= total_dry_mass * 1e9
    i_counter += 1

i_counter = 0
for [time, filename, key] in time_filename_list2:
    print time, filename, key
    ncf = Scientific.IO.NetCDF.NetCDFFile(filename)
    particles = partmc.aero_particle_array_t(ncf)
    env_state = partmc.env_state_t(ncf)
    ncf.close()

    total_number = sum(1/particles.comp_vols)
    total_dry_mass = sum(particles.mass(exclude = ["H2O"])/particles.comp_vols)

    array_num[i_counter,1]= total_number
    array_mass[i_counter,1]= total_dry_mass * 1e9
    i_counter += 1

i_counter = 0
for [time, filename, key] in time_filename_list3:
    print time, filename, key
    ncf = Scientific.IO.NetCDF.NetCDFFile(filename)
    particles = partmc.aero_particle_array_t(ncf)
    env_state = partmc.env_state_t(ncf)
    ncf.close()

    total_number = sum(1/particles.comp_vols)
    total_dry_mass = sum(particles.mass(exclude = ["H2O"])/particles.comp_vols)

    array_num[i_counter,2]= total_number
    array_mass[i_counter,2]= total_dry_mass * 1e9
    i_counter += 1

print array_num.shape, array_mass.shape, time_array.shape
num_avg = np.average(array_num, axis = 1)
mass_avg = np.average(array_mass, axis = 1)

plt.clf()
plt.plot(time_array[:], array_num[:,0], 'r')
plt.plot(time_array[:], array_num[:,1], 'b')
plt.plot(time_array[:], array_num[:,2], 'g')
plt.plot(time_array[:], num_avg[:], 'k')
plt.xlabel("time ")
plt.ylabel("number concentration in m^{-3}")
plt.title("weighted -2")
fig = plt.gcf()
fig.savefig("figs/number_coarse_w-2.pdf")

plt.clf()
plt.plot(time_array[:], array_mass[:,0], 'r')
plt.plot(time_array[:], array_mass[:,1], 'b')
plt.plot(time_array[:], array_mass[:,2], 'g')
plt.plot(time_array[:], mass_avg[:], 'k')
plt.xlabel("time ")
plt.ylabel("mass concentration in \mu g m^{-3}")
plt.title("weighted -2")
fig = plt.gcf()
fig.savefig("figs/mass_coarse_w-2.pdf")




    



