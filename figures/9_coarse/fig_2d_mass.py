#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc

def make_plot(dir_name,in_filename,out_filename):
    ncf = Scientific.IO.NetCDF.NetCDFFile(dir_name+in_filename)
    particles = pmc_data_nc.aero_particle_array_t(ncf)
    ncf.close()

    bc = particles.mass(include = ["NO3"])
    dry_mass = particles.mass(exclude = ["H2O"])
    bc_frac = bc / dry_mass

    dry_diameter = particles.dry_diameter()

    x_axis = pmc_data_nc.pmc_log_axis(min=1e-9,max=1e-5,n_bin=70)
    y_axis = pmc_data_nc.pmc_linear_axis(min=0,max=1.,n_bin=50)

    hist2d = pmc_data_nc.histogram_2d(dry_diameter, bc_frac, x_axis, y_axis, weights = particles.mass(include = ["NO3"])/particles.comp_vol)
    plt.clf()
    plt.pcolor(x_axis.edges(), y_axis.edges(), hist2d.transpose(),norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("linear")
    plt.axis([x_axis.min, x_axis.max, y_axis.min, y_axis.max])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("NO3 mass fraction")
    plt.title("NO3 mass distribution")
    cbar = plt.colorbar()
    cbar.set_label("number density (m^{-3})")
    fig = plt.gcf()
    fig.savefig(out_filename)

dir_name = "../../scenarios/5_coarse/out/"

#filename_in = "urban_plume_nc_0001_00000001.nc"
#filename_out = "figs/2d_nc_bc_01.pdf"
#make_plot(dir_name, filename_in, filename_out)

#filename_in = "urban_plume_nc_0001_00000012.nc"
#filename_out = "figs/2d_nc_bc_12-w2.pdf"
#make_plot(dir_name, filename_in, filename_out)

#filename_in = "urban_plume_nc_0001_00000018.nc"
#filename_out = "figs/2d_nc_bc_18-w2.pdf"
#make_plot(dir_name, filename_in, filename_out)

#filename_in = "urban_plume_nc_0001_00000024.nc"
#filename_out = "figs/2d_nc_bc_24.pdf"
#make_plot(dir_name, filename_in, filename_out)

#filename_in = "urban_plume_wc_0001_00000001.nc"
#filename_out = "figs/2d_wc_bc_01.pdf"
#make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000012.nc"
filename_out = "figs/2d_wc_no3_12-w2.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000018.nc"
filename_out = "figs/2d_wc_no3_18-w2.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000024.nc"
filename_out = "figs/2d_wc_no3_24-w2.pdf"
make_plot(dir_name, filename_in, filename_out)

