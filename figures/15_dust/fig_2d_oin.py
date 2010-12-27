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

def make_plot(in_filename,out_filename,title):
    ncf = scipy.io.netcdf.netcdf_file(in_filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    ncf.close()

    oin = particles.masses(include = ["OIN"])
    dry_mass = particles.masses(exclude = ["H2O"])
    oin_frac = oin / dry_mass

    dry_diameters = particles.dry_diameters()

    x_axis = partmc.log_grid(min=1e-8,max=1e-4,n_bin=100)
    y_axis = partmc.linear_grid(min=0,max=1,n_bin=50)

    hist2d = partmc.histogram_2d(dry_diameters, oin_frac, x_axis, y_axis, weights = 1/particles.comp_vols)
    plt.clf()
    plt.pcolor(x_axis.edges(), y_axis.edges(), hist2d.transpose(),norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("linear")
    plt.axis([1e-8, 1e-4, 0, 1])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("OIN mass fraction")
    plt.grid(True)
    plt.clim(1e6, 1e10)
    cbar = plt.colorbar()
    cbar.set_label(r"number density ($\rm m^{-3}$)")
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

for hour in range(1, 25):
    print "hour = ", hour
    
    filename_in1 = "../../scenarios/7_dust/out/urban_plume_wc_0001_000000%02d.nc" % hour
    filename_out1 = "figs/2d_oin_100K_wei-2_emit_%02d.pdf" % (hour-1)
    titel = "%02d hours" % (hour-1)
    print filename_in1
    print filename_out1
    print titel

    make_plot(filename_in1, filename_out1, titel)
#    make_plot(filename_in2, filename_out2, titel)
#    make_plot(filename_in3, filename_out3, titel)
#    make_plot(filename_in4, filename_out4, titel)


