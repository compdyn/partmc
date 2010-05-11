#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

def make_plot(dir_name,in_filename,out_filename):
    ncf = scipy.io.netcdf.netcdf_file(dir_name+in_filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    ncf.close()

    bc = particles.masses(include = ["BC"])
    dry_mass = particles.masses(exclude = ["H2O"])
    bc_frac = bc / dry_mass

    dry_diameters = particles.dry_diameters()

    x_axis = partmc.log_grid(min=1e-10,max=1e-6,n_bin=70)
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
    fig = plt.gcf()
    fig.savefig(out_filename)

dir_name = "../../scenarios/4_nucleate/out/"

filename_in = "urban_plume_wc_0001_00000001.nc"
filename_out = "figs/2d_wc_bc_01.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000037.nc"
filename_out = "figs/2d_wc_bc_37.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000073.nc"
filename_out = "figs/2d_wc_bc_73.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000109.nc"
filename_out = "figs/2d_wc_bc_109.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000145.nc"
filename_out = "figs/2d_wc_bc_145.pdf"
make_plot(dir_name, filename_in, filename_out)


