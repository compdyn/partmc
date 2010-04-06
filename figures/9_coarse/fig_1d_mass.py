#!/usr/bin/env python

import Scientific.IO.NetCDF
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
        ncf = Scientific.IO.NetCDF.NetCDFFile(in_dir+file)
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
#    plt.axis([1e-10, 1e-4, 1e4, 1e11])
    plt.xlabel("dry diameter (m)")
    plt.ylabel(" mass concentration (kg m^{-3})")
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

dir_name = "../../scenarios/5_weighted/out/"
for hour in range(1, 26):
    print "hour = ", hour
    files = []
    print 'files ', files
    for i_loop in range (0, config.i_loop_max):
        filename_in = "urban_plume_wc_100K_wei-3_00%02d_000000%02d.nc" % ((i_loop+1), hour)
        print i_loop, filename_in
        files.append(filename_in)
    print files
    filename_out = "figs/1d_mass_100K_wei-3_%02d.pdf" % hour
    title = '100K wei-3, %02d hour' % (hour-1)
    print filename_out
    make_plot(dir_name, files, title, filename_out)
