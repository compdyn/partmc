#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc
import mpl_helper

def make_plot(in_dir, in_filename1, in_filename2, in_filename3, out_filename):
    print in_filename1, in_filename2, in_filename3
    ncf = scipy.io.netcdf.netcdf_file(in_dir+in_filename1, 'r')
    particles1 = partmc.aero_particle_array_t(ncf)
    ncf.close()
    ncf = scipy.io.netcdf.netcdf_file(in_dir+in_filename2, 'r')
    particles2 = partmc.aero_particle_array_t(ncf)
    ncf.close()
    ncf = scipy.io.netcdf.netcdf_file(in_dir+in_filename3, 'r')
    particles3 = partmc.aero_particle_array_t(ncf)
    ncf.close()

    x_axis = partmc.log_grid(min=1e-10,max=1e-4,n_bin=50)
    x_centers = x_axis.centers() 

    dry_diameters1 = particles1.dry_diameters()
    dry_diameters2 = particles2.dry_diameters()
    dry_diameters3 = particles3.dry_diameters()

    hist1 = partmc.histogram_1d(dry_diameters1, x_axis, weights = particles1.masses(exclude=["H2O"])/ particles1.comp_vols)
    hist2 = partmc.histogram_1d(dry_diameters2, x_axis, weights = particles2.masses(exclude=["H2O"]) / particles2.comp_vols)
    hist3 = partmc.histogram_1d(dry_diameters3, x_axis, weights = particles3.masses(exclude=["H2O"]) / particles3.comp_vols)
 
    plt.clf()
    plt.loglog(x_axis.centers(), hist1, label = 'initial')
    plt.loglog(x_axis.centers(), hist2, label = '6 hours')
    plt.loglog(x_axis.centers(), hist3, label = '12 hours') 
    plt.legend(loc = 'center right')
    plt.axis([5e-9, 1e-4, 1e-14, 1e-7])
    plt.grid(True)
    plt.xlabel("dry diameter (m)")
    plt.ylabel(r"mass concentration ($\rm kg \, m^{-3}$)")
    fig = plt.gcf()
    fig.savefig(out_filename)

    
#dir_name = "../../scenarios/7_dust/out_wei-3_emit/"
dir_name = "../../scenarios/7_dust/out/"

filename_in1 = "urban_plume_wc_0001_00000001.nc" 
filename_in2 = "urban_plume_wc_0001_00000007.nc" 
filename_in3 = "urban_plume_wc_0001_00000013.nc" 

filename_out = "figs/1d_wc_mass_100K_wei-2_emit.pdf" 

make_plot(dir_name, filename_in1, filename_in2, filename_in3, filename_out)

