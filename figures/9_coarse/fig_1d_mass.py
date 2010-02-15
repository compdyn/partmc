#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

#def make_plot(in_dir, in_filename1, in_filename2, in_filename3, out_filename):
def make_plot(in_dir, in_filename1, out_filename):
    print in_filename1
    ncf1 = Scientific.IO.NetCDF.NetCDFFile(in_dir+in_filename1)
    particles1 = partmc.aero_particle_array_t(ncf1)
    ncf1.close()

#    ncf2 = Scientific.IO.NetCDF.NetCDFFile(in_dir+in_filename2)
#    particles2 = partmc.aero_particle_array_t(ncf2)
#    ncf2.close()

#    ncf3 = Scientific.IO.NetCDF.NetCDFFile(in_dir+in_filename3)
#    particles3 = partmc.aero_particle_array_t(ncf3)
#    ncf3.close()

    x_axis = partmc.log_grid(min=1e-10,max=1e-4,n_bin=100)
    x_centers = x_axis.centers() 

    dry_diameters1 = particles1.dry_diameters()
#    dry_diameters2 = particles2.dry_diameters()
#    dry_diameters3 = particles3.dry_diameters()

    hist1 = partmc.histogram_1d(dry_diameters1, x_axis, weights = particles1.masses() / particles1.comp_vols)
#    hist2 = partmc.histogram_1d(dry_diameters2, x_axis, weights = particles2.masses() / particles2.comp_vols)
#    hist3 = partmc.histogram_1d(dry_diameters3, x_axis, weights = particles3.masses() / particles3.comp_vols)

    plt.clf()
    plt.loglog(x_axis.centers(), hist1, 'r')
#    plt.loglog(x_axis.centers(), hist2, 'b')
#    plt.loglog(x_axis.centers(), hist3, 'g')
#    plt.axis([1e-10, 1e-4, 1e7, 1e15])
    plt.xlabel("dry diameter (m)")
    plt.ylabel("mass density (kg m^{-3})")
    plt.title("100K flat")
    fig = plt.gcf()
    fig.savefig(out_filename)

dir_name = "../../scenarios/5_weighted/out/"
for hour in range(1, 26):
    print "hour = ", hour
    filename_in1 = "urban_plume_wc_0001_000000%02d.nc" % hour
#    filename_in2 = "urban_plume_wc_0002_000000%02d.nc" % hour
#    filename_in3 = "urban_plume_wc_0003_000000%02d.nc" % hour
    filename_out = "figs/1d_wc_mass_%02d.pdf" % hour
#    make_plot(dir_name, filename_in1, filename_in2, filename_in3, filename_out)
    make_plot(dir_name, filename_in1, filename_out)
