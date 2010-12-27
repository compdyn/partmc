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

#netcdf_dir = "../../scenarios/7_dust/out/"
netcdf_dir = "../../scenarios/2_urban_plume2/out/"
netcdf_pattern = "urban_plume_wc_0001_(.*).nc"
time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)

num_array = np.zeros([len(time_filename_list),3])
i_counter = 0

for [time, filename, key] in time_filename_list:
    print time, filename, key
    ncf = scipy.io.netcdf.netcdf_file(filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    env_state = partmc.env_state_t(ncf)
    ncf.close()

    num_total = sum(1 / particles.comp_vols)
    dry_diameters = particles.dry_diameters()
    is_small = (dry_diameters < 5e-8)

    num_small = sum(1/particles.comp_vols[is_small])

    num_array[i_counter,0] = time
    num_array[i_counter,1] = num_total / 1e6
    num_array[i_counter,2] = num_small / 1e6

    i_counter += 1

print "maximum number conc ", max(num_array[:,1])
plt.plot(num_array[:,0]/3600., num_array[:,1], label = 'total number')
plt.plot(num_array[:,0]/3600., num_array[:,2], label = 'number small')
plt.xlim(0, 24)
plt.xticks([0, 6, 12, 18, 24])
plt.grid(True)
plt.xlabel("time in hours")
plt.ylabel(r"number concentration in $\rm cm^{-3}$")
fig = plt.gcf()
#fig.savefig("figs/number_time_100K_wei-2_emit.pdf")
fig.savefig("figs/number_time_up2_compare.pdf")


    



