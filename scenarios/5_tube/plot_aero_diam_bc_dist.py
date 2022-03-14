#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import matplotlib
import scipy.io, numpy

for (filename, index) in mpl_helper.get_filename_list('out_as2suc/', r'urban_plume_([0-9]+)_process\.nc'):
    (figure, axes, cbar_axes) = mpl_helper.make_fig(left_margin=0.7, right_margin=1, colorbar=True)

    ncf = scipy.io.netcdf_file(filename)
    diam_edges = ncf.variables["diam_edges"].data.copy() * 1e6
    oc_frac_edges = ncf.variables["oc_frac_edges"].data.copy() * 100
    diam_oc_dist = ncf.variables["diam_oc_dist"].data.copy() * 1e-6

    p = axes.pcolor(diam_edges, oc_frac_edges, diam_oc_dist,
                    norm = matplotlib.colors.LogNorm(vmin=1e3, vmax=1e5), linewidths = 0.1)

    axes.set_xscale("log")
    axes.set_xlabel(r"dry diameter $D_{\rm dry}$ / $\rm \mu m$")
    axes.set_xlim(1e-2, 1e0)

    axes.set_yscale("linear")
    axes.set_ylabel(r"OC mass fraction $w_{\rm OC}$ / \%")
    axes.set_ylim(0, 100)

    axes.grid(True)
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                           orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number conc. $n(D_{\rm dry},w_{\rm OC})$ / $\rm cm^{-3}$")

    out_filename = "out_as2suc/urban_plume_diam_oc_dist_%s.pdf" % index
    print("Writing %s" % out_filename)
    figure.savefig(out_filename)
    matplotlib.pyplot.close(figure)
