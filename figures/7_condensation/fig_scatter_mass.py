#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
import subprocess
import os

matplotlib.use("PDF")
matplotlib.use('Agg')

import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc

def make_plot(in_dir, in_filename_wc, in_filename_nc, title, out_filename_wc, out_filename_nc):
    print 'file ', in_dir+in_filename_wc
    
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_dir+in_filename_wc)
    particles_wc = pmc_data_nc.aero_particle_array_t(ncf)
    ncf.close()
    
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_dir+in_filename_nc)
    particles_nc = pmc_data_nc.aero_particle_array_t(ncf)
    ncf.close()

#    so4_wc =  particles_wc.mass(include = ["SO4"])/particles_wc.aero_data.molec_weight[0]
#    nh4_wc =  particles_wc.mass(include = ["NH4"])/particles_wc.aero_data.molec_weight[3] 
#    no3_wc =  particles_wc.mass(include = ["NO3"])/particles_wc.aero_data.molec_weight[1]

    so4_nc =  particles_nc.mass(include = ["SO4"])/particles_nc.aero_data.molec_weight[0]
    nh4_nc =  particles_nc.mass(include = ["NH4"])/particles_nc.aero_data.molec_weight[3] 
    no3_nc =  particles_nc.mass(include = ["NO3"])/particles_nc.aero_data.molec_weight[1]  

#    plt.scatter(nh4_wc, 2*so4_wc+no3_wc)
    
#    a = plt.gca() # gets the axis
#    a.set_xscale("log") # x axis log
#    a.set_yscale("log") # y axis log
#    plt.axis([1e-25, 1e-15, 1e-25, 1e-15]) # axis limit

#    plt.title(title)
#    fig = plt.gcf()
#    fig.savefig(out_filename_wc)

    plt.scatter(nh4_nc, 2*so4_nc+no3_nc)
    
    a = plt.gca() # gets the axis
    a.set_xscale("log") # x axis log
    a.set_yscale("log") # y axis log
    plt.axis([1e-25, 1e-15, 1e-25, 1e-15]) # axis limit

    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename_nc)

for counter in range(10,11):
    print "counter = ", counter
    
    in_dir = "../../urban_plume2/out_no_nh3/"
    in_filename_wc = "urban_plume_wc_0001_000000%02d.nc" % counter
    in_filename_nc = "urban_plume_nc_0001_000000%02d.nc" % counter
    title = " %02d hours" % (counter-1)
    out_filename_wc = "figs/scatter_mass_wc_%02d.pdf" % counter
    out_filename_nc = "figs/scatter_mass_nc_%02d.pdf" % counter

    print title

    make_plot(in_dir, in_filename_wc, in_filename_nc, title, out_filename_wc, out_filename_nc)

