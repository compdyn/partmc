#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc
const = pmc_data_nc.load_constants("../../src/constants.f90")

netcdf_dir = "../../urban_plume_nucleate/out/"
netcdf_pattern = "urban_plume_nc_0001_(.*).nc"
time_filename_list = pmc_data_nc.get_time_filename_list(netcdf_dir, netcdf_pattern)

size_dist_array = np.zeros([len(time_filename_list),100])
times = np.zeros([len(time_filename_list)])
i_counter = 0
diam_axis = pmc_data_nc.pmc_log_axis(min=1e-10,max=1e-6,n_bin=100)
diam_axis_edges = diam_axis.edges()

for [time, filename, key] in time_filename_list:
    print time, filename, key
    ncf = Scientific.IO.NetCDF.NetCDFFile(filename)
    particles = pmc_data_nc.aero_particle_array_t(ncf)
    ncf.close()

    dry_diameter = particles.dry_diameter()
    hist = pmc_data_nc.histogram_1d(dry_diameter, diam_axis, weights = 1 / particles.comp_vol)
    size_dist_array[i_counter,:] = hist
    times[i_counter] = time
    i_counter += 1

np.savetxt("data/banana_size_dist.txt", size_dist_array)
np.savetxt("data/banana_diam.txt", diam_axis_edges)
np.savetxt("data/banana_times.txt", times)
    



