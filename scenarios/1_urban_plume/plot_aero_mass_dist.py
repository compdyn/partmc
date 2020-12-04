#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import matplotlib
import scipy.io, numpy

for (filename, index) in mpl_helper.get_filename_list('out/', r'urban_plume_ch_([0-9]+)_process\.nc'):
    (figure, axes) = mpl_helper.make_fig(left_margin=0.7, right_margin=1)

    ncf = scipy.io.netcdf_file(filename)
    diam_edges = ncf.variables["diam_edges"].data.copy() * 1e6
    diam_bc_dist = ncf.variables["mass_bc_dist"].data.copy() * 1e9
    diam_so4_dist = ncf.variables["mass_so4_dist"].data.copy()* 1e9
    diam_no3_dist = ncf.variables["mass_no3_dist"].data.copy()* 1e9
    diam_nh4_dist = ncf.variables["mass_nh4_dist"].data.copy()* 1e9
    
    axes.plot(diam_edges[:-1], diam_bc_dist, label="BC")
    axes.plot(diam_edges[:-1], diam_so4_dist, label="SO4")
    axes.plot(diam_edges[:-1], diam_no3_dist, label="NO3")
    axes.plot(diam_edges[:-1], diam_nh4_dist, label="NH4")
    
    axes.set_xscale("log")
    axes.set_xlabel(r"dry diameter $D_{\rm dry}$ / $\rm \mu m$")
    axes.set_xlim(1e-1, 1e1)

    axes.set_yscale("log")
    axes.set_ylabel(r"dm/dlogDp / $\rm \mu$g/$\rm m^{-3}$")
    #axes.set_ylim(1e-10, 1e-4)

    axes.grid(True)
    axes.legend()
    out_filename = "out/urban_plume_bc_mass_dist_%s.pdf" % index
    print("Writing %s" % out_filename)
    figure.savefig(out_filename)
    matplotlib.pyplot.close(figure)
