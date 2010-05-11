#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
import subprocess 
import os

matplotlib.use("PDF")
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

def make_plot(in_filename,out_filename):
    print in_filename
    ncf = scipy.io.netcdf.netcdf_file(in_filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    env_state = partmc.env_state_t(ncf)
    ncf.close()

    dry_diameters = particles.dry_diameters()*1e6
    s_crit = (particles.critical_rel_humids(env_state) - 1)*100
    x_axis = partmc.log_grid(min=1e-2,max=1e0,n_bin=70)
    y_axis = partmc.log_grid(min=1e-3,max=1e2,n_bin=50)

    hist2d = partmc.histogram_2d(dry_diameters, s_crit, x_axis, y_axis, weights = 1/particles.comp_vols)
    plt.clf()
    plt.pcolor(x_axis.edges(), y_axis.edges(), hist2d.transpose(),norm = matplotlib.colors.LogNorm(), linewidths = 0.1)

    title = partmc.time_of_day_string(env_state)
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("log")
    plt.axis([x_axis.min, x_axis.max, y_axis.min, y_axis.max])
    plt.xlabel("dry diameter (\mu m)")
    plt.ylabel("critical supersaturation (%)")
    cbar = plt.colorbar()
    plt.clim(1e8, 1e10)
    cbar.set_label("number density (m^{-3})")
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

for hour in range(1, 1442):
    print "hour = ", hour

    filename_in1 = "../../scenarios/2_urban_plume2/out/urban_plume_wc_0001_0000%04d.nc" % hour
    filename_out1 = "figs/2d_scrit_ref_%04d.png" % (hour-1)
    print filename_in1
    print filename_out1

    make_plot(filename_in1, filename_out1)
 
