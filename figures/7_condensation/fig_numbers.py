#!/usr/bin/env python

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc
const = pmc_data_nc.load_constants("../../src/constants.f90")

in_dir = "../../new_cond/out/"
out_filename = "figs/numbers_wc.pdf" 

time_array = np.linspace(0,47,48)
frac_array = np.zeros((4,len(time_array))) # run x time
run_list = ["ref", "comp", "size", "both"]

for k in range(0,4):
    run = run_list[k]
    for counter in range(0,48):
        in_filename_start = "cond_wc_%02d_%s_0001_00000001.nc" % (counter+1, run)
        in_filename_end = "cond_wc_%02d_%s_0001_00000601.nc" % (counter+1, run)
        ncf = Scientific.IO.NetCDF.NetCDFFile(in_dir+in_filename_start)
        particles_start = pmc_data_nc.aero_particle_array_t(ncf)
        ncf.close()
        ncf = Scientific.IO.NetCDF.NetCDFFile(in_dir+in_filename_end)
        particles_end = pmc_data_nc.aero_particle_array_t(ncf)
        ncf.close()

        num_start = sum(1/particles_start.comp_vol)
        num_end = sum(1/particles_end.comp_vol)

        frac_array[k,counter] = np.double(num_end) / np.double(num_start)
        print particles_start.comp_vol, particles_end.comp_vol
        print 'num ', run, counter, num_end, num_start, frac_array[k,counter]
    print run, np.max(frac_array[k,:]), np.min(frac_array[k,:])

plt.clf()
plt.plot(time_array, frac_array[0,:], 'r-', label = 'ref')
plt.plot(time_array, frac_array[1,:], 'r.', label = 'comp')
plt.plot(time_array, frac_array[2,:], 'r--', label = 'size')
plt.plot(time_array, frac_array[3,:], 'r-.', label = 'both')
plt.legend(loc = 'lower right')
plt.xlabel("time ")
plt.ylabel("number_end/number_start")
fig = plt.gcf()
fig.savefig(out_filename)

