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

def make_plot(in_dir, in_files, title, out_filename, error):
    x_axis = partmc.log_grid(min=1e-8,max=1e-5,n_bin=3)
    x_centers = x_axis.centers()
    counter = 0
    hist_array = np.zeros([len(x_centers),config.i_loop_max])
    error = np.zeros([3])
    for file in in_files:
        ncf = scipy.io.netcdf.netcdf_file(in_dir+file, 'r')
        particles = partmc.aero_particle_array_t(ncf)
        ncf.close() 

        dry_diameters = particles.dry_diameters()
        hist = partmc.histogram_1d(dry_diameters, x_axis, weights = 1 / particles.comp_vols)
        hist_array[:,counter] = hist
        counter = counter+1
    plt.clf()
    for i_loop in range(0,config.i_loop_max):
        plt.loglog(x_axis.centers(), hist_array[:,i_loop], 'k')
        plt.errorbar(x_axis.centers(), np.average(hist_array,axis = 1), np.std(hist_array, axis = 1))
        avg = np.average(hist_array, axis = 1)
        std = np.std(hist_array, axis = 1)
        error = std/avg
	print 'avg and std ', avg, std, error
    plt.axis([1e-8, 1e-5, 1e4, 1e11])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("number density (m^{-3})")
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

dir_name = "../../scenarios/5_weighted/out/"
error = np.zeros([3])
array_error = np.zeros([3,25])

for hour in range(1, 26):
    print "hour = ", hour
    files = []
    for i_loop in range (0, config.i_loop_max):
        filename_in = "urban_plume_wc_10K_wei-3_00%02d_000000%02d.nc" % ((i_loop+1), hour)
        print i_loop, filename_in
        files.append(filename_in)
    filename_out = "figs/1d_10K_wei-3_%02d_3bin.pdf" % hour
    title = '10K wei-3, %02d hour' % (hour-1)
    make_plot(dir_name, files, title, filename_out, error)
    array_error[:,hour-1] = error
    print 'error ', error 
print 'array_error ', array_error
