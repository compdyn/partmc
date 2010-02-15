#!/usr/bin/env python

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

netcdf_dir = "../../scenarios/4_nucleate/out/"
netcdf_pattern = "urban_plume_wc_0001_(.*).nc"

time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)

size_dist_array = np.zeros([len(time_filename_list),100])
times = np.zeros([len(time_filename_list)])
i_counter = 0
diam_axis = partmc.log_grid(min=1e-10,max=1e-6,n_bin=100)
diam_axis_edges = diam_axis.edges()

for [time, filename, key] in time_filename_list:
    print time, filename, key
    ncf = Scientific.IO.NetCDF.NetCDFFile(filename)
    particles = partmc.aero_particle_array_t(ncf)
    ncf.close()

    dry_diameters = particles.dry_diameters()
    hist = partmc.histogram_1d(dry_diameters, diam_axis, weights = 1 / particles.comp_vols)
    size_dist_array[i_counter,:] = hist
    times[i_counter] = time
    i_counter += 1

np.savetxt("data/banana_size_dist_wc.txt", size_dist_array)
np.savetxt("data/banana_diam_wc.txt", diam_axis_edges)
np.savetxt("data/banana_times_wc.txt", times)
    



