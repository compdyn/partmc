#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

def make_plot(in_dir, in_filename, out_filename):
    print in_filename
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_dir+in_filename)
    particles = partmc.aero_particle_array_t(ncf)
    ncf.close()

    x_axis = partmc.log_grid(min=1e-10,max=1e-4,n_bin=100)
    x_centers = x_axis.centers() 

    dry_diameters = particles.dry_diameters()

    hist = partmc.histogram_1d(dry_diameters, x_axis, weights = particles.masses() / particles.comp_vols)

    plt.clf()
    plt.loglog(x_axis.centers(), hist)
#    plt.axis([1e-10, 1e-4, 1e7, 1e15])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("mass density (kg m^{-3})")
    fig = plt.gcf()
    fig.savefig(out_filename)

dir_name = "../../scenarios/4_nucleate/out/"

#filename_in = "urban_plume_wc_0001_00000001.nc"
#filename_out = "figs/1d_wc_mass_001.pdf"
#make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000005.nc"
filename_out = "figs/1d_wc_mass_005_w.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000013.nc"
filename_out = "figs/1d_wc_mass_013_w.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000016.nc"
filename_out = "figs/1d_wc_mass_016_w.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000018.nc"
filename_out = "figs/1d_wc_mass_018_w.pdf"
make_plot(dir_name, filename_in, filename_out)

filename_in = "urban_plume_wc_0001_00000024.nc"
filename_out = "figs/1d_wc_mass_024_w.pdf"
make_plot(dir_name, filename_in, filename_out)
