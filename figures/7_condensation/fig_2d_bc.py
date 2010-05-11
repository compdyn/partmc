#!/usr/bin/env python2.5

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
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

    dry_diameters = particles.dry_diameters()

    x_axis = partmc.log_grid(min=1e-8,max=1e-6,n_bin=70)
    y_axis = partmc.linear_grid(min=0,max=0.8,n_bin=40)

    hist2d = partmc.histogram_2d(dry_diameters, bc_frac, x_axis, y_axis, weights = 1/particles.comp_vols)
    plt.clf()
    plt.pcolor(x_axis.edges(), y_axis.edges(), hist2d.transpose(),norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("linear")
    plt.axis([x_axis.min, x_axis.max, y_axis.min, y_axis.max])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("BC mass fraction")
    cbar = plt.colorbar()
    cbar.set_label("number density (m^{-3})")
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

for hour in range(1, 25):
    print "hour = ", hour
    
    filename_in1 = "../../scenarios/1_urban_plume/out/urban_plume_wc_0001_000000%02d.nc" % hour
    filename_in2 = "../../scenarios/3_condense/start/urban_plume_comp_wc_0001_000000%02d.nc" % hour
    filename_in3 = "../../scenarios/3_condense/start/urban_plume_size_wc_0001_000000%02d.nc" % hour
    filename_in4 = "../../scenarios/3_condense/start/urban_plume_both_wc_0001_000000%02d.nc" % hour
    filename_out1 = "figs/2d_bc_urban_plume_%02d.pdf" % (hour-1)
    filename_out2 = "figs/2d_bc_comp_%02d.pdf" % (hour-1)
    filename_out3 = "figs/2d_bc_size_%02d.pdf" % (hour-1)
    filename_out4 = "figs/2d_bc_both_%02d.pdf" % (hour-1)
    titel = "%02d hours" % (hour-1)
    print filename_in1
    print filename_out1
    print titel

    make_plot(filename_in1, filename_out1, titel)
#    make_plot(filename_in2, filename_out2, titel)
#    make_plot(filename_in3, filename_out3, titel)
#    make_plot(filename_in4, filename_out4, titel)


