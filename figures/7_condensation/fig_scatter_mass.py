#!/usr/bin/env python

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

def make_plot(in_dir, in_filename_wc, in_filename_nc, out_filename, title):
    print in_dir
    
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_filename_wc)
    particles_wc = pmc_data_nc.aero_particle_array_t(ncf)
    ncf.close()
    
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_filename_nc)
    particles_nc = pmc_data_nc.aero_particle_array_t(ncf)
    ncf.close()

    so4_wc =  particles_wc.mass(include = ["SO4"])
    nh4_wc =  particles_wc.mass(include = ["NH4"]) 
    scatter(so4_wc,nh4_wc) # colorscale log
    a = gca() # gets the axis
    a.set_xscale("log") # x axis log
    a.set_yscale("log") # y axis log
    axis([1e-25, 1e-15, 1e-25, 1e-15]) # axis limit

    plt.title(title)
    fig = plt.gcf()
    fig.savefig(out_filename)

for counter in range(1,49):
    print "counter = ", counter
    
    dir_name = "../../urban_plume2/out"
    filename_in1 = "cond_%02d_ref_0001_.*.nc" % counter
    filename_in2 = "cond_%02d_comp_0001_.*.nc" % counter
    filename_in3 = "cond_%02d_size_0001_.*.nc" % counter
    filename_in4 = "cond_%02d_both_0001_.*.nc" % counter

    filename_out1 = "figs/env_ref_%02d.png" % (counter-1)
    filename_out2 = "figs/env_comp_%02d.png" % (counter-1)
    filename_out3 = "figs/env_size_%02d.png" % (counter-1)
    filename_out4 = "figs/env_both_%02d.png" % (counter-1)

    title = " %02d hours" % (counter-1)

    print dir_name, title
    print filename_in1, filename_out1
    print filename_in2, filename_out2
    print filename_in3, filename_out3
    print filename_in4, filename_out4

    make_plot(dir_name,filename_in1, filename_out1, title, 0, counter-1)
    make_plot(dir_name,filename_in2, filename_out2, title, 1, counter-1)
    make_plot(dir_name,filename_in3, filename_out3, title, 2, counter-1)
    make_plot(dir_name,filename_in4, filename_out4, title, 3, counter-1)

np.savetxt("data/maximum_ss.txt", maximum_ss)
