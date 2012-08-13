#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

def make_plot(in_filename,out_filename):
    ncf = scipy.io.netcdf.netcdf_file(in_filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    ncf.close()

    bc = particles.masses(include = ["BC"])
    dry_mass = particles.masses(exclude = ["H2O"])
    bc_frac = bc / dry_mass

    wet_diameters = particles.diameters()

    x_axis = partmc.log_grid(min=1e-8,max=1e-4,n_bin=100)
    y_axis = partmc.linear_grid(min=0,max=0.8,n_bin=40)

    hist2d = partmc.histogram_2d(wet_diameters, bc_frac, x_axis, y_axis, weights = particles.num_concs)
    plt.clf()
    plt.pcolor(x_axis.edges(), y_axis.edges(), hist2d.transpose(),norm = matplotlib.colors.LogNorm(vmin=1e8,vmax=1e12), linewidths = 0.1)
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("linear")
    plt.axis([x_axis.min, x_axis.max, y_axis.min, y_axis.max])
    plt.xlabel("wet diameter (m)")
    plt.ylabel("BC mass fraction")
    cbar = plt.colorbar()
    cbar.set_label("number density (m^{-3})")
    fig = plt.gcf()
    fig.savefig(out_filename)

filename_in1 = "../../scenarios/8_condense_p/out_00000004/condense_0001_00000001.nc" 
filename_in2 = "../../scenarios/8_condense_p/out_00000004/condense_0001_00000031.nc" 
filename_in3 = "../../scenarios/8_condense_p/out_00000004/condense_0001_00000061.nc" 
filename_out1 = "figs/2d_bc_wet_diam_00000004_01.pdf" 
filename_out2 = "figs/2d_bc_wet_diam_00000004_31.pdf" 
filename_out3 = "figs/2d_bc_wet_diam_00000004_61.pdf" 

make_plot(filename_in1, filename_out1)
make_plot(filename_in2, filename_out2)
make_plot(filename_in3, filename_out3)

filename_in1 = "../../scenarios/8_condense_p/out_00000008/condense_0001_00000001.nc" 
filename_in2 = "../../scenarios/8_condense_p/out_00000008/condense_0001_00000031.nc" 
filename_in3 = "../../scenarios/8_condense_p/out_00000008/condense_0001_00000061.nc" 
filename_out1 = "figs/2d_bc_wet_diam_00000008_01.pdf" 
filename_out2 = "figs/2d_bc_wet_diam_00000008_31.pdf" 
filename_out3 = "figs/2d_bc_wet_diam_00000008_61.pdf" 

make_plot(filename_in1, filename_out1)
make_plot(filename_in2, filename_out2)
make_plot(filename_in3, filename_out3)

filename_in1 = "../../scenarios/8_condense_p/out_00000013/condense_0001_00000001.nc" 
filename_in2 = "../../scenarios/8_condense_p/out_00000013/condense_0001_00000031.nc" 
filename_in3 = "../../scenarios/8_condense_p/out_00000013/condense_0001_00000061.nc" 
filename_out1 = "figs/2d_bc_wet_diam_00000013_01.pdf" 
filename_out2 = "figs/2d_bc_wet_diam_00000013_31.pdf" 
filename_out3 = "figs/2d_bc_wet_diam_00000013_61.pdf" 

make_plot(filename_in1, filename_out1)
make_plot(filename_in2, filename_out2)
make_plot(filename_in3, filename_out3)

filename_in1 = "../../scenarios/8_condense_p/out_00000020/condense_0001_00000001.nc" 
filename_in2 = "../../scenarios/8_condense_p/out_00000020/condense_0001_00000031.nc" 
filename_in3 = "../../scenarios/8_condense_p/out_00000020/condense_0001_00000061.nc" 
filename_out1 = "figs/2d_bc_wet_diam_00000020_01.pdf" 
filename_out2 = "figs/2d_bc_wet_diam_00000020_31.pdf" 
filename_out3 = "figs/2d_bc_wet_diam_00000020_61.pdf" 

make_plot(filename_in1, filename_out1)
make_plot(filename_in2, filename_out2)
make_plot(filename_in3, filename_out3)
