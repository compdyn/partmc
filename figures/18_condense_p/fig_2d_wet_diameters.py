#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
sys.path.append("../../tool")
import mpl_helper
import matplotlib
matplotlib.use("PDF")
#import matplotlib.pyplot as plt
import partmc

def make_plot(in_filename,out_filename):
    (figure, axes, cbar_axes) = mpl_helper.make_fig(left_margin=0.7, right_margin=1, colorbar=True)
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
    p =  axes.pcolor(x_axis.edges(), y_axis.edges(), hist2d.transpose(),norm = matplotlib.colors.LogNorm(vmin=1e8,vmax=1e12), linewidths = 0.1)

    axes.set_xscale("log")
    axes.set_xlabel(r"dry diameter $D_{\rm dry}$ / $\rm \mu m$")
    axes.set_xlim(1e-8, 1e-4)

    axes.set_yscale("linear")
    axes.set_ylabel(r"BC mass fraction $w_{\rm BC}$ / \%")
    axes.set_ylim(0, 0.8)

    axes.grid(True)
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                           orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number conc. $n(D_{\rm dry},w_{\rm BC})$ / $\rm cm^{-3}$")

    figure.savefig(out_filename)
    print out_filename


filename_in1 = "../../scenarios/8_condense_p/out_we/out_04/condense_0001_00000001.nc" 
filename_in2 = "../../scenarios/8_condense_p/out_we/out_04/condense_0001_00000031.nc" 
filename_in3 = "../../scenarios/8_condense_p/out_we/out_04/condense_0001_00000061.nc" 
filename_out1 = "figs/2d_bc_wet_diam_we_04_01.pdf" 
filename_out2 = "figs/2d_bc_wet_diam_we_04_31.pdf" 
filename_out3 = "figs/2d_bc_wet_diam_we_04_61.pdf" 

make_plot(filename_in1, filename_out1)
make_plot(filename_in2, filename_out2)
make_plot(filename_in3, filename_out3)

filename_in1 = "../../scenarios/8_condense_p/out_we/out_08/condense_0001_00000001.nc" 
filename_in2 = "../../scenarios/8_condense_p/out_we/out_08/condense_0001_00000031.nc" 
filename_in3 = "../../scenarios/8_condense_p/out_we/out_08/condense_0001_00000061.nc" 
filename_out1 = "figs/2d_bc_wet_diam_we_08_01.pdf" 
filename_out2 = "figs/2d_bc_wet_diam_we_08_31.pdf" 
filename_out3 = "figs/2d_bc_wet_diam_we_08_61.pdf" 

make_plot(filename_in1, filename_out1)
make_plot(filename_in2, filename_out2)
make_plot(filename_in3, filename_out3)

filename_in1 = "../../scenarios/8_condense_p/out_we/out_13/condense_0001_00000001.nc" 
filename_in2 = "../../scenarios/8_condense_p/out_we/out_13/condense_0001_00000031.nc" 
filename_in3 = "../../scenarios/8_condense_p/out_we/out_13/condense_0001_00000061.nc" 
filename_out1 = "figs/2d_bc_wet_diam_we_13_01.pdf" 
filename_out2 = "figs/2d_bc_wet_diam_we_13_31.pdf" 
filename_out3 = "figs/2d_bc_wet_diam_we_13_61.pdf" 

make_plot(filename_in1, filename_out1)
make_plot(filename_in2, filename_out2)
make_plot(filename_in3, filename_out3)

filename_in1 = "../../scenarios/8_condense_p/out_we/out_20/condense_0001_00000001.nc" 
filename_in2 = "../../scenarios/8_condense_p/out_we/out_20/condense_0001_00000031.nc" 
filename_in3 = "../../scenarios/8_condense_p/out_we/out_20/condense_0001_00000061.nc" 
filename_out1 = "figs/2d_bc_wet_diam_we_20_01.pdf" 
filename_out2 = "figs/2d_bc_wet_diam_we_20_31.pdf" 
filename_out3 = "figs/2d_bc_wet_diam_we_20_61.pdf" 

make_plot(filename_in1, filename_out1)
make_plot(filename_in2, filename_out2)
make_plot(filename_in3, filename_out3)
