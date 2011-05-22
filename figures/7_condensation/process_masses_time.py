#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc
import mpl_helper

import config

netcdf_dir = "../../scenarios/1_urban_plume/out/"
netcdf_pattern = "urban_plume_nc_0001_(.*).nc"
time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)

mass_array = np.zeros([len(time_filename_list),5])
i_counter = 0

for [time, filename, key] in time_filename_list:
    print time, filename, key
    ncf = scipy.io.netcdf.netcdf_file(filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    env_state = partmc.env_state_t(ncf)
    ncf.close()

    bc = particles.masses(include = ["BC"])
    bc_total = sum(particles.masses(include = ["BC"])/particles.comp_vols)

    is_bc = (bc > 0)
    is_not_bc = (bc == 0)

    inorg = particles.masses(include = ["SO4", "NO3", "NH4"])
    inorg_on_bc = inorg[is_bc]
    inorg_not_on_bc = inorg[is_not_bc]

    print "inorg ", inorg, inorg_on_bc
    print "sums ", sum(inorg/particles.comp_vols), sum(inorg_on_bc/particles.comp_vols[is_bc])

    mass_array[i_counter,0] = time
    mass_array[i_counter,1] = sum(inorg/particles.comp_vols) * 1e9
    mass_array[i_counter,2] = sum(inorg_on_bc/particles.comp_vols[is_bc]) * 1e9
    mass_array[i_counter,3] = sum(inorg_not_on_bc/particles.comp_vols[is_not_bc]) * 1e9

    mass_array[i_counter,4] = bc_total * 1e9

    i_counter += 1

plt.plot(mass_array[:,0]/3600., mass_array[:,1], label = 'total inorg')
plt.plot(mass_array[:,0]/3600., mass_array[:,2], label = 'total inorg on BC')
plt.plot(mass_array[:,0]/3600., mass_array[:,3], label = 'total inorg on non-BC')
plt.plot(mass_array[:,0]/3600., mass_array[:,4], label = 'total BC')

plt.legend(loc = 'center right')
plt.title("Inorganic mass and BC mass (no coag)")
plt.grid(True)
plt.xlabel("time in hours")
plt.ylabel(r"mass in $\rm \mu g \, m^{-3}$")

fig = plt.gcf()
fig.savefig("figs/bc_inorg_time_nc.pdf")



    



