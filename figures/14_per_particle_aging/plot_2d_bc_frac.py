#!/usr/bin/env python

import os, sys
import scipy.io
import numpy as np

sys.path.append("/Users/nriemer/subversion/partmc/trunk/tool")
import partmc
import mpl_helper
import matplotlib
import pickle

def make_plot(dir_name,in_filename,out_filename):
    ncf = scipy.io.netcdf.netcdf_file(dir_name+in_filename, 'r')
    particles = partmc.aero_particle_array_t(ncf)
    ncf.close()

    bc = particles.masses(include = ["BC"])
    dry_mass = particles.masses(exclude = ["H2O"])
    bc_frac = bc / dry_mass

    dry_diameters = particles.dry_diameters() * 1e6

    x_axis = partmc.log_grid(min=1e-3,max=1e1,n_bin=100)
    y_axis = partmc.linear_grid(min=0,max=0.8,n_bin=40)

    hist2d = partmc.histogram_2d(dry_diameters, bc_frac, x_axis, y_axis, weights = 1/particles.comp_vols)

    hist2d = hist2d * 1e-6
    print hist2d[36,:]
    (figure, axes_array, cbar_axes_array) = mpl_helper.make_fig_array(1,1, figure_width=5,
									  top_margin=0.5, bottom_margin=0.45,
									  left_margin=0.65, right_margin=1,
									  vert_sep=0.3, horiz_sep=0.3,
									  colorbar="shared", colorbar_location="right")

    axes = axes_array[0][0]
    cbar_axes = cbar_axes_array[0]
    p = axes.pcolor(x_axis.edges(), y_axis.edges(), hist2d.transpose(), 
		    norm = matplotlib.colors.LogNorm(vmin=1e3, vmax=1e5), linewidths = 0.1)
    axes.set_xscale("log")
    axes.set_yscale("linear")
    axes.set_ylabel(r"BC mass fraction $w_{\rm BC}$")
    axes.set_xlabel(r"dry diameter $D$/ $\rm \mu m$")
    axes.set_ylim(0,0.8)
    axes.set_xlim(5e-3, 1e0)
    axes.grid(True)
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
			   orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number conc. $n(D,w_{\rm BC})$ / $\rm cm^{-3}$")
    
    mpl_helper.remove_fig_array_axes(axes_array)
    figure.savefig(out_filename)

dir_name = "../../scenarios/1_urban_plume/out/"

filename_in = "urban_plume_wc_0001_00000120.nc"
filename_out = "figs/2d_wc_bc_120.pdf"
make_plot(dir_name, filename_in, filename_out)


			
