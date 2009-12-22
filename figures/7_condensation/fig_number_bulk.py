#!/usr/bin/env python

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc
import config
const = pmc_data_nc.load_constants("../../src/constants.f90")

netcdf_dir = "../../urban_plume2/out/"
netcdf_pattern = "urban_plume_wc_0001_(.*).nc"
time_filename_list = pmc_data_nc.get_time_filename_list(netcdf_dir, netcdf_pattern)

array_wc = np.zeros([len(time_filename_list),4])
array_nc = np.zeros([len(time_filename_list),4])

i_counter = 0
for [time, filename, key] in time_filename_list:
    print time, filename, key
    ncf = Scientific.IO.NetCDF.NetCDFFile(filename)
    particles = pmc_data_nc.aero_particle_array_t(ncf)
    env_state = pmc_data_nc.env_state_t(ncf)
    ncf.close()

    total_number = sum(1/particles.comp_vol)
    total_dry_mass = sum(particles.mass(exclude = ["H2O"])/particles.comp_vol)
    total_mass = sum(particles.mass()/particles.comp_vol)
    array_wc[i_counter,0]= time / 3600.
    array_wc[i_counter,1]= total_number
    array_wc[i_counter,2]= total_dry_mass
    array_wc[i_counter,3]= total_mass
    i_counter += 1

netcdf_dir = "../../urban_plume2/out/"
netcdf_pattern = "urban_plume_nc_0001_(.*).nc"
time_filename_list = pmc_data_nc.get_time_filename_list(netcdf_dir, netcdf_pattern)

i_counter = 0
for [time, filename, key] in time_filename_list:
    print time, filename, key
    ncf = Scientific.IO.NetCDF.NetCDFFile(filename)
    particles = pmc_data_nc.aero_particle_array_t(ncf)
    env_state = pmc_data_nc.env_state_t(ncf)
    ncf.close()

    total_number = sum(1/particles.comp_vol)
    total_dry_mass = sum(particles.mass(exclude = ["H2O"])/particles.comp_vol)
    total_mass = sum(particles.mass()/particles.comp_vol)

    array_nc[i_counter,0]= time / 3600.
    array_nc[i_counter,1]= total_number
    array_nc[i_counter,2]= total_dry_mass
    array_nc[i_counter,3]= total_mass
    i_counter += 1

plt.clf()
plt.plot(array_wc[:,0], array_wc[:,1], 'r-', label = 'wc')
plt.plot(array_nc[:,0], array_nc[:,1], 'r--', label = 'nc')
plt.legend(loc = 'upper right')
plt.xlabel("time ")
plt.ylabel("number concentration in m^{-3}")
fig = plt.gcf()
fig.savefig("figs/number.pdf")

plt.clf()
plt.plot(array_wc[:,0], array_wc[:,2], 'r-', label = 'dry wc')
plt.plot(array_wc[:,0], array_wc[:,3], 'b-', label = 'total wc')
plt.plot(array_nc[:,0], array_nc[:,2], 'r--', label = 'dry nc')
plt.plot(array_nc[:,0], array_nc[:,3], 'b--', label = 'total nc')
plt.legend(loc = 'upper right')
plt.xlabel("time ")
plt.ylabel("mass concentration in kg m^{-3}")
fig = plt.gcf()
fig.savefig("figs/mass.pdf")

    



