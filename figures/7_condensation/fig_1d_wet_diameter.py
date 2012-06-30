#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

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

    x_axis = partmc.log_grid(min=1e-10,max=1e-4,n_bin=30)
    x_centers = x_axis.centers() 

    wet_diameters1 = particles1.diameters()
    wet_diameters2 = particles2.diameters()
    wet_diameters3 = particles3.diameters()

    hist1 = partmc.histogram_1d(wet_diameters1, x_axis, weights = particles1.num_concs)
    hist2 = partmc.histogram_1d(wet_diameters2, x_axis, weights = particles2.num_concs)
    hist3 = partmc.histogram_1d(wet_diameters3, x_axis, weights = particles3.num_concs)
 
    plt.clf()
    plt.semilogx(x_axis.centers(), hist1, label = '0 min')
    plt.semilogx(x_axis.centers(), hist2, label = '5 mins')
    plt.semilogx(x_axis.centers(), hist3, label = '10 mins') 
    plt.ylim(0, 1.4e10)
    plt.legend(loc = 'upper left')
    plt.xlabel("wet diameter (m)")
    plt.ylabel("number density (m^{-3})")
    plt.grid(True)
    fig = plt.gcf()
    fig.savefig(out_filename)

dir_name = "../../scenarios/8_condense_p/out/"

filename_in1 = "condense_0001_00000001.nc" 
filename_in2 = "condense_0001_00000031.nc" 
filename_in3 = "condense_0001_00000061.nc" 

filename_out = "figs/1d_wet_diameter_ne.pdf" 
print dir_name
print filename_in1, filename_in2, filename_in3
print filename_out

make_plot(dir_name, filename_in1, filename_in2, filename_in3, filename_out)

