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
    print in_filename
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

make_plot("../../new_cond/start_up2/urban_plume_wc_0001_00000031.nc","figs/2d_bc_30.pdf","30 hours")
make_plot("../../new_cond/start_up2/urban_plume_wc_0001_00000037.nc","figs/2d_bc_36.pdf","36 hours")
make_plot("../../new_cond/start_up2/urban_plume_wc_0001_00000043.nc","figs/2d_bc_42.pdf","42 hours")
make_plot("../../new_cond/start_up2/urban_plume_wc_0001_00000049.nc","figs/2d_bc_48.pdf","48 hours")

make_plot("../../new_cond/start_up2_comp/urban_plume_wc_0001_00000031.nc","figs/2d_bc_comp_30.pdf","30 hours")
make_plot("../../new_cond/start_up2_comp/urban_plume_wc_0001_00000037.nc","figs/2d_bc_comp_36.pdf","36 hours")
make_plot("../../new_cond/start_up2_comp/urban_plume_wc_0001_00000043.nc","figs/2d_bc_comp_42.pdf","42 hours")
make_plot("../../new_cond/start_up2_comp/urban_plume_wc_0001_00000049.nc","figs/2d_bc_comp_48.pdf","48 hours")

make_plot("../../new_cond/start_up2_size/urban_plume_wc_0001_00000031.nc","figs/2d_bc_size_30.pdf","30 hours")
make_plot("../../new_cond/start_up2_size/urban_plume_wc_0001_00000037.nc","figs/2d_bc_size_36.pdf","36 hours")
make_plot("../../new_cond/start_up2_size/urban_plume_wc_0001_00000043.nc","figs/2d_bc_size_42.pdf","42 hours")
make_plot("../../new_cond/start_up2_size/urban_plume_wc_0001_00000049.nc","figs/2d_bc_size_48.pdf","48 hours")

make_plot("../../new_cond/start_up2_both/urban_plume_wc_0001_00000031.nc","figs/2d_bc_both_30.pdf","30 hours")
make_plot("../../new_cond/start_up2_both/urban_plume_wc_0001_00000037.nc","figs/2d_bc_both_36.pdf","36 hours")
make_plot("../../new_cond/start_up2_both/urban_plume_wc_0001_00000043.nc","figs/2d_bc_both_42.pdf","42 hours")
make_plot("../../new_cond/start_up2_both/urban_plume_wc_0001_00000049.nc","figs/2d_bc_both_48.pdf","48 hours")
