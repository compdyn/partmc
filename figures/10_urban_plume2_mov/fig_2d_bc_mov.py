#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
matplotlib.use('Agg')

import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

def make_plot(in_filename,out_filename,title):
    ncf = scipy.io.netcdf.netcdf_file(in_filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    ncf.close()

    bc = particles.masses(include = ["BC"])
    dry_mass = particles.masses(exclude = ["H2O"])
    bc_frac = bc / dry_mass

    wet_diameters = particles.diameters()
    dry_diameters = particles.dry_diameters() * 1e6

    x_axis = partmc.log_grid(min=1e-2,max=1e0,n_bin=90)
    y_axis = partmc.linear_grid(min=0,max=0.8,n_bin=40)

    hist2d = partmc.histogram_2d(dry_diameters, bc_frac, x_axis, y_axis, weights = 1/particles.comp_vols)
    plt.clf()
    plt.pcolor(x_axis.edges(), y_axis.edges(), hist2d.transpose(),norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("linear")
    plt.axis([x_axis.min, x_axis.max, y_axis.min, y_axis.max])
    plt.xlabel("dry diameter (\mu m)")
    plt.ylabel("BC dry mass fraction")
    cbar = plt.colorbar()
    plt.clim(1e8, 1e11)
    cbar.set_label("number density (m^{-3})")
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

for counter in range(1, 1442):
    print "counter = ",  counter
    
    filename_in1 = "../../scenarios/2_urban_plume2/out/urban_plume_wc_0001_0000%04d.nc" % counter
    filename_out1 = "figs/2d_bc_%04d.png" % (counter-1)
    titel = "%02d minutes" % (counter-1)
    print filename_in1
    print filename_out1
    print titel

    make_plot(filename_in1, filename_out1, titel)
#    make_plot(filename_in2, filename_out2, titel)
#    make_plot(filename_in3, filename_out3, titel)
#    make_plot(filename_in4, filename_out4, titel)

