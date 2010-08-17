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

    x_axis = partmc.log_grid(min=1e-8,max=1e-4,n_bin=90)
    y_axis = partmc.linear_grid(min=0,max=0.8,n_bin=40)

    hist2d = partmc.histogram_2d(wet_diameters, bc_frac, x_axis, y_axis, weights = 1/particles.comp_vols)
    plt.clf()
    plt.pcolor(x_axis.edges(), y_axis.edges(), hist2d.transpose(),norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("linear")
    plt.axis([x_axis.min, x_axis.max, y_axis.min, y_axis.max])
    plt.xlabel("wet diameter (m)")
    plt.ylabel("BC mass fraction")
    cbar = plt.colorbar()
    plt.clim(1e9, 1e12)
    cbar.set_label("number density (m^{-3})")
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

for counter in range(1, 602):
    print "counter = ",  counter
    
    filename_in1 = "/home/ching1/subversion/partmc/trunk/scenarios/3_condense/out/cond_02_ref_0001_00000%03d.nc" % counter
    filename_out1 = "figs/2d_bc_ref_%03d.png" % (counter-1)
    titel = "%02d seconds" % (counter-1)
    print filename_in1
    print filename_out1
    print titel

    make_plot(filename_in1, filename_out1, titel)
#    make_plot(filename_in2, filename_out2, titel)
#    make_plot(filename_in3, filename_out3, titel)
#    make_plot(filename_in4, filename_out4, titel)

