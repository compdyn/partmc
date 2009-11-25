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

make_plot("../../urban_plume2/out/urban_plume_nc_0001_00000002.nc","fig1.pdf","1 hours")
make_plot("../../urban_plume2/out/urban_plume_nc_0001_00000013.nc","fig12.pdf","12 hours")
make_plot("../../urban_plume2/out/urban_plume_nc_0001_00000025.nc","fig24.pdf","24 hours")
make_plot("../../urban_plume2/out/urban_plume_nc_0001_00000049.nc","fig48.pdf","48 hours")

