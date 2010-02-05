#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc

def make_plot(in_filename,out_filename,title):
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_filename)
    particles = pmc_data_nc.aero_particle_array_t(ncf)
    ncf.close()

    bc = particles.mass(include = ["BC"])
    dry_mass = particles.mass(exclude = ["H2O"])
    bc_frac = bc / dry_mass

    dry_diameter = particles.dry_diameter()

    x_axis = pmc_data_nc.pmc_log_axis(min=1e-8,max=1e-6,n_bin=70)
    y_axis = pmc_data_nc.pmc_linear_axis(min=0,max=0.8,n_bin=40)

    hist2d = pmc_data_nc.histogram_2d(dry_diameter, bc_frac, x_axis, y_axis, weights = 1/particles.comp_vol)
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
    
    filename_in1 = "../../urban_plume/out/urban_plume_wc_0001_000000%02d.nc" % hour
    filename_in2 = "../../new_cond/start/urban_plume_comp_wc_0001_000000%02d.nc" % hour
    filename_in3 = "../../new_cond/start/urban_plume_size_wc_0001_000000%02d.nc" % hour
    filename_in4 = "../../new_cond/start/urban_plume_both_wc_0001_000000%02d.nc" % hour
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


