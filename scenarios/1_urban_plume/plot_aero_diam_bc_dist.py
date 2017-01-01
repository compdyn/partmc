#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import matplotlib
import partmc
import scipy.io, numpy

for (filename, index) in partmc.get_filename_list('out/', r'urban_plume_([0-9]+)_process\.nc'):
    (figure, axes, cbar_axes) = mpl_helper.make_fig(left_margin=0.7, right_margin=1, colorbar=True)

    ncf = scipy.io.netcdf_file(filename)
    diam_edges = ncf.variables["diam_edges"].data.copy() * 1e6
    bc_frac_edges = ncf.variables["bc_frac_edges"].data.copy() * 100
    diam_bc_dist = ncf.variables["diam_bc_dist"].data.copy() * 1e-6

    p = axes.pcolor(diam_edges, bc_frac_edges, diam_bc_dist,
                    norm = matplotlib.colors.LogNorm(vmin=1e3, vmax=1e5), linewidths = 0.1)

    axes.set_xscale("log")
    axes.set_xlabel(r"dry diameter $D_{\rm dry}$ / $\rm \mu m$")
    axes.set_xlim(1e-2, 1e0)

    axes.set_yscale("linear")
    axes.set_ylabel(r"BC mass fraction $w_{\rm BC}$ / \%")
    axes.set_ylim(0, 80)

    axes.grid(True)
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                           orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number conc. $n(D_{\rm dry},w_{\rm BC})$ / $\rm cm^{-3}$")

    out_filename = "out/urban_plume_diam_bc_dist_%s.pdf" % index
    print("Writing %s" % out_filename)
    figure.savefig(out_filename)
    matplotlib.pyplot.close(figure)
