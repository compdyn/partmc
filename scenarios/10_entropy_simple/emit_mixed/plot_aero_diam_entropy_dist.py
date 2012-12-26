#!/usr/bin/env python

import sys, os
sys.path.append("../../../tool")
import mpl_helper
import matplotlib
import partmc
import scipy.io, numpy

for (filename, index) in partmc.get_filename_list('out/', r'urban_plume_([0-9]+)_process\.nc'):
    (figure, axes, cbar_axes) = mpl_helper.make_fig(left_margin=0.7, right_margin=1, colorbar=True)

    ncf = scipy.io.netcdf_file(filename)
    diam_edges = ncf.variables["diam_edges"].data * 1e6
    entropy_edges = ncf.variables["entropy_edges"].data
    diam_entropy_dist = ncf.variables["diam_entropy_dist"].data * 1e-6

    p = axes.pcolor(diam_edges, entropy_edges, diam_entropy_dist,
                    norm = matplotlib.colors.LogNorm(vmin=1e3, vmax=1e5), linewidths = 0.1)

    axes.set_xscale("log")
    axes.set_xlabel(r"dry diameter $D_{\rm dry}$ / $\rm \mu m$")
    axes.set_xlim(1e-2, 1e0)

    axes.set_yscale("linear")
    axes.set_ylabel(r"entropy")
    axes.set_ylim(0, 2)

    axes.grid(True)
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                           orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number conc. $n(D_{\rm dry},s)$ / $\rm cm^{-3}$")

    out_filename = "out/urban_plume_diam_entropy_distdry_%s.pdf" % index
    figure.savefig(out_filename)
    print out_filename
