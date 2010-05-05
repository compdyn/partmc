#!/usr/bin/env python

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
matplotlib.use('Agg')

import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

def make_plot(in_filename,out_filename):
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_filename)
    particles = partmc.aero_particle_array_t(ncf)
    env_state = partmc.env_state_t(ncf)
    ncf.close()

    bc = particles.masses(include = ["BC"])
    dry_mass = particles.masses(exclude = ["H2O"])
    bc_frac = bc / dry_mass

    so4 = particles.masses(include = ["SO4"])
    inorg = particles.masses(include = ["SO4", "NO3", "NH4"])
    inorg_frac = inorg / dry_mass 

    kappas = particles.kappas()
 
    wet_diameters = particles.diameters()
    dry_diameters = particles.dry_diameters() * 1e6

    x_axis = partmc.log_grid(min=1e-2,max=1e0,n_bin=90)
    y_axis = partmc.linear_grid(min=0,max=0.8,n_bin=40)

    vals = partmc.multival_2d(dry_diameters, bc_frac, kappas, x_axis, y_axis, rand_arrange=False)

    vals_pos = np.ma.masked_less_equal(vals, 0)
    vals_zero = np.ma.masked_not_equal(vals, 0)
    
    plt.clf()
    if vals_zero.count() > 0:
       plt.pcolor(x_axis.edges(), y_axis.edges(), vals_zero.transpose(), cmap=matplotlib.cm.gray, linewidths = 0.1)
    
    if vals_pos.count() > 0:
       plt.pcolor(x_axis.edges(), y_axis.edges(), vals_pos.transpose(), linewidths = 0.1)

    title = partmc.time_of_day_string(env_state)
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("linear")
    plt.axis([x_axis.min, x_axis.max, y_axis.min, y_axis.max])
    plt.xlabel("dry diameter (\mu m)")
    plt.ylabel("BC dry mass fraction")
    cbar = plt.colorbar()
    plt.clim(0, 0.6)
    cbar.set_label("kappa")
    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

for counter in range(1, 1442):
    print "counter = ",  counter
    
    filename_in1 = "../../scenarios/2_urban_plume2/out/urban_plume_wc_0001_0000%04d.nc" % counter
    filename_out1 = "figs/2d_kappa_%04d.png" % (counter-1)
    print filename_in1
    print filename_out1

    make_plot(filename_in1, filename_out1)

