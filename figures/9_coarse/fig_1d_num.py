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

def make_plot(in_dir, in_files, title, out_filename):
    x_axis = partmc.log_grid(min=1e-10,max=1e-4,n_bin=100)
    x_centers = x_axis.centers()
    counter = 0
    hist_array = np.zeros([len(x_centers),config.i_loop_max])
    for file in in_files:
        ncf = scipy.io.netcdf.netcdf_file(in_dir+file, 'r')
        particles = partmc.aero_particle_array_t(ncf)
        ncf.close() 

        dry_diameters = particles.dry_diameters()
        hist = partmc.histogram_1d(dry_diameters, x_axis, weights = 1 / particles.comp_vols)
        hist_array[:,counter] = hist
        counter = counter+1
#    hist_array_av = np.average(hist_array,axis = 1)
#    hist_array_std = np.std(hist_array, axis = 1)
#    hist_array_std_clipped = np.minimum(hist_array_std, hist_array_av - 1)
#    e_bars = np.vstack((hist_array_std_clipped, hist_array_std))

    hist_array_gav = np.exp(np.average(np.log(hist_array),axis = 1))
    hist_array_gstd = np.exp(np.std(np.log(hist_array), axis = 1))
    e_bar_top = hist_array_gav * hist_array_gstd 
    e_bar_bottom = hist_array_gav / hist_array_gstd
    e_bars = np.vstack((hist_array_gav - e_bar_bottom, e_bar_top - hist_array_gav))


    plt.clf()
#    for i_loop in range(0,config.i_loop_max):
#        plt.loglog(x_axis.centers(), hist_array[:,i_loop], 'k')
    a = plt.gca() # gets the axis
    a.set_xscale("log") # x axis log
    a.set_yscale("log") # y axis log
    plt.errorbar(x_axis.centers(), hist_array_gav, e_bars)
    plt.axis([5e-9, 5e-6, 1e4, 1e11])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("number density (m^{-3})")
#    plt.title(title)
    plt.grid(True)
    fig = plt.gcf()
    fig.savefig(out_filename)

dir_name = "../../scenarios/5_weighted/out/"
for case in ["10K_wei+1", "10K_flat", "10K_wei-1", "10K_wei-2", "10K_wei-3", "10K_wei-4" ]:
#for case in ["100K_wei+1", "100K_flat", "100K_wei-1", "100K_wei-2", "100K_wei-3", "100K_wei-4" ]:
    for hour in range(12, 13):
        print "hour = ", hour
        files = []
        print 'files ', files
        for i_loop in range (0, config.i_loop_max):
            filename_in = "urban_plume_wc_%s_0%03d_000000%02d.nc" % (case, (i_loop+1), hour)
            print i_loop, filename_in
            files.append(filename_in)
        print files
        filename_out = "figs/100loop/1d_%s_%02d.pdf" % (case, hour)
        title = '%s, %02d hour' % (case, hour-1)
        print filename_out
        make_plot(dir_name, files, title, filename_out)
