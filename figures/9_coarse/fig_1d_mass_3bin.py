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

def make_plot(in_dir, in_files, title, out_filename, hour):
    x_axis = partmc.log_grid(min=1e-8,max=1e-5,n_bin=3)
    x_centers = x_axis.centers()
    counter = 0
    hist_array = np.zeros([len(x_centers),config.i_loop_max])
    for file in in_files:
        ncf = scipy.io.netcdf.netcdf_file(in_dir+file, 'r')
        particles = partmc.aero_particle_array_t(ncf)
        ncf.close() 

        dry_diameters = particles.dry_diameters()
        hist = partmc.histogram_1d(dry_diameters, x_axis, weights = particles.masses(exclude=["H2O"]) / particles.comp_vols)
        hist_array[:,counter] = hist
        counter = counter+1
    plt.clf()
    for i_loop in range(0,config.i_loop_max):
        plt.loglog(x_axis.centers(), hist_array[:,i_loop], 'k')
        plt.errorbar(x_axis.centers(), np.average(hist_array,axis = 1), np.std(hist_array, axis = 1))
        avg = np.average(hist_array, axis = 1)
        std = np.std(hist_array, axis = 1)
        error[:,hour-1] = std/avg
    print 'error ', hour, error
    plt.xlabel("dry diameter (m)")
    plt.ylabel("mass concentration (kg m^{-3})")
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

dir_name = "../../scenarios/5_weighted/out/"
error = np.zeros([3,25])
avg_error = np.zeros([5,3])
i_counter = 0
for counter in ["flat", "wei-1", "wei-2", "wei-3", "wei-4"]:
    print 'counter', counter
    for hour in range(1, 26):
        print "hour = ", hour
        files = []
        for i_loop in range (0, config.i_loop_max):
            filename_in = "urban_plume_wc_10K_%s_00%02d_000000%02d.nc" % (counter, (i_loop+1), hour)
            files.append(filename_in)
        filename_out = "figs/1d_mass_10K_%s_%02d_3bin.pdf" % (counter, hour)
        title = '10K %s, %02d hour' % (counter, hour-1)
        make_plot(dir_name, files, title, filename_out, hour)
    avg_error[i_counter,:] = np.average(error, axis = 1)
    i_counter = i_counter + 1
f1 = 'data/mass_1d_3bin_10K.txt'
np.savetxt(f1, avg_error)

