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

netcdf_dir = "../../scenarios/7_dust/out_wei-3_emit/"
netcdf_pattern = "urban_plume_wc_0001_(.*).nc"
time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)

mass_array = np.zeros([len(time_filename_list),6])
num_array = np.zeros([len(time_filename_list)])
i_counter = 0

for [time, filename, key] in time_filename_list:
    print time, filename, key
    ncf = scipy.io.netcdf.netcdf_file(filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    env_state = partmc.env_state_t(ncf)
    ncf.close()

    num_total = sum(1 / particles.comp_vols)
    bc = particles.masses(include = ["BC"])
    bc_total = sum(particles.masses(include = ["BC"])/particles.comp_vols)

    oin = particles.masses(include = ["OIN"])
    oin_total = sum(particles.masses(include = ["OIN"])/particles.comp_vols)

    is_oin = (oin > 0)
    is_not_oin = (oin == 0)

    inorg = particles.masses(include = ["SO4", "NO3", "NH4"])
    inorg_on_oin = inorg[is_oin]
    inorg_not_on_oin = inorg[is_not_oin]

    mass_array[i_counter,0] = time
    mass_array[i_counter,1] = sum(inorg/particles.comp_vols) * 1e9
    mass_array[i_counter,2] = sum(inorg_on_oin/particles.comp_vols[is_oin]) * 1e9
    mass_array[i_counter,3] = sum(inorg_not_on_oin/particles.comp_vols[is_not_oin]) * 1e9

    mass_array[i_counter,4] = bc_total * 1e9
    mass_array[i_counter,5] = oin_total * 1e9

    i_counter += 1

plt.plot(mass_array[:,0]/3600., mass_array[:,1], label = 'total inorg')
plt.plot(mass_array[:,0]/3600., mass_array[:,2], label = 'total inorg on OIN')
plt.plot(mass_array[:,0]/3600., mass_array[:,3], label = 'total inorg on non-OIN')
plt.plot(mass_array[:,0]/3600., mass_array[:,4], label = 'total BC')
plt.plot(mass_array[:,0]/3600., mass_array[:,5], label = 'total OIN')
plt.xlim(0, 24)
plt.xticks([0, 6, 12, 18, 24])

plt.legend(loc = 'upper right')
#plt.title("Inorganic mass and OIN mass")
plt.grid(True)
plt.xlabel("time in hours")
plt.ylabel(r"mass concentration in $\rm \mu g \, m^{-3}$")

fig = plt.gcf()
fig.savefig("figs/bc_inorg_time_10K_wei-3_emit.pdf")




    



