#!/usr/bin/env python

import sys, os
sys.path.append("../../../tool")
import mpl_helper
import matplotlib
import partmc
import scipy.io, numpy

for (filename, index) in partmc.get_filename_list('out/', r'urban_plume_0001_([0-9]+)_process\.nc'):
    (figure, axes, cbar_axes) = mpl_helper.make_fig(left_margin=0.7, right_margin=1, colorbar=True)

    ncf = scipy.io.netcdf_file(filename)
    entropies = ncf.variables["entropies"].data
    dry_diameters = ncf.variables["dry_diam"].data * 1e6
    scs = ncf.variables["sc"].data * 100
    bc_fracs = ncf.variables["bc_frac"].data * 100
    h2o_fracs = ncf.variables["h2o_frac"].data * 100
    oc_fracs = ncf.variables["oc_frac"].data * 100
    no3_fracs = ncf.variables["no3_frac"].data * 100
    so4_fracs = ncf.variables["so4_frac"].data * 100
    
    least_ct = ncf.variables["least_ct"].data
    greatest_ct = ncf.variables["greatest_ct"].data
    
    p = axes.scatter(dry_diameters,entropies,c=greatest_ct)#,norm=matplotlib.colors.LinNorm()) 

    axes.set_xscale("log")
    axes.set_xlabel(r"dry diameter $D_{\rm dry}$ / $\rm \mu m$")
    axes.set_xlim(1e-2, 1e0)

    axes.set_yscale("linear")
    axes.set_ylabel(r"entropy")
    axes.set_ylim(0, 2)

    axes.grid(True)

    cbar = figure.colorbar(p, cax=cbar_axes, orientation='vertical') 
#format=matplotlib.ticker.LogFormatterMathtext(),
                           
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"greatest create times")

    out_filename = "out/urban_plume_diam_entropy_scatter_greatest_ct_%s.pdf" % index
    figure.savefig(out_filename)
    print out_filename
