#!/usr/bin/env python2.5

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

def make_plot(in_filename,out_filename,time,title):
    ncf = scipy.io.netcdf.netcdf_file(in_filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    ncf.close()

    age = abs(particles.least_create_times / 3600. - time)
    dry_diameters = particles.dry_diameters()

    x_axis = partmc.log_grid(min=1e-8,max=1e-6,n_bin=70)
    y_axis = partmc.linear_grid(min=0, max = 48, n_bin=49)

    hist2d = partmc.histogram_2d(dry_diameters, age, x_axis, y_axis, weights = 1/particles.comp_vols)

    plt.clf()
    plt.pcolor(x_axis.edges(), y_axis.edges(), hist2d.transpose(),norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("linear")
    plt.axis([x_axis.min, x_axis.max, y_axis.min, y_axis.max])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("age (h)")
    cbar = plt.colorbar()
    cbar.set_label("number density (m^{-3})")
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

for hour in range(49, 50):
    print "hour = ", hour
    time = hour - 1
    filename_in1 = "../../scenarios/2_urban_plume2/out/urban_plume_nc_0001_000000%02d.nc" % hour
    filename_out1 = "figs/2d_age_nc_%02d.pdf" % (hour-1)
    titel = "%02d hours" % (hour-1)
    print filename_in1
    print filename_out1
    print titel

    make_plot(filename_in1, filename_out1, time, titel)


