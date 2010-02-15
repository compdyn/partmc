#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

def make_plot(dir_name,in_filename,out_filename):
    ncf = Scientific.IO.NetCDF.NetCDFFile(dir_name+in_filename)
    particles = partmc.aero_particle_array_t(ncf)
    ncf.close()

    bc = particles.masses(include = ["H2O"])
    dry_mass = particles.masses(exclude = ["H2O"])
    bc_frac = bc / dry_mass

    dry_diameters = particles.dry_diameters()

    x_axis = partmc.log_grid(min=1e-9,max=1e-5,n_bin=70)
    y_axis = partmc.linear_grid(min=0,max=1.,n_bin=50)

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

dir_name = "../../scenarios/5_weighted/out/"

filename_in = "urban_plume_nc_0001_00000001.nc"
filename_out = "figs/2d_nc_h2o_01.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_nc_0001_00000012.nc"
filename_out = "figs/2d_nc_h2o_12.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_nc_0001_00000018.nc"
filename_out = "figs/2d_nc_h2o_18.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_nc_0001_00000024.nc"
filename_out = "figs/2d_nc_h2o_24.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000001.nc"
filename_out = "figs/2d_wc_h2o_01.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000012.nc"
filename_out = "figs/2d_wc_h2o_12.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000018.nc"
filename_out = "figs/2d_wc_h2o_18.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000024.nc"
filename_out = "figs/2d_wc_h2o_24.pdf"
make_plot(dir_name, filename_in, filename_out)

