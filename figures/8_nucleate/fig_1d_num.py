#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc

def make_plot(in_dir, in_filename, out_filename):
    print in_filename
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_dir+in_filename)
    particles = pmc_data_nc.aero_particle_array_t(ncf)
    ncf.close()

    x_axis = pmc_data_nc.pmc_log_axis(min=1e-10,max=1e-4,n_bin=100)
    x_centers = x_axis.centers() 

    dry_diameter = particles.dry_diameter()

    hist = pmc_data_nc.histogram_1d(dry_diameter, x_axis, weights = 1 / particles.comp_vol)

    plt.clf()
    plt.loglog(x_axis.centers(), hist)
    plt.axis([1e-10, 1e-4, 1e7, 1e15])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("number density (m^{-3})")
    fig = plt.gcf()
    fig.savefig(out_filename)

dir_name = "../../urban_plume/out/"

filename_in = "urban_plume_wc_0001_00000001.nc"
filename_out = "figs/1d_wc_num_001.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000005.nc"
filename_out = "figs/1d_wc_num_005.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000013.nc"
filename_out = "figs/1d_wc_num_013.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000016.nc"
filename_out = "figs/1d_wc_num_016.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000018.nc"
filename_out = "figs/1d_wc_num_018.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000024.nc"
filename_out = "figs/1d_wc_num_024.pdf"
make_plot(dir_name, filename_in, filename_out)
