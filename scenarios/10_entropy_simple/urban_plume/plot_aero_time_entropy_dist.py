#!/usr/bin/env python

import sys, os
sys.path.append("../../../tool")
import mpl_helper
import matplotlib
import partmc
import scipy.io, numpy

filename = 'out/urban_plume2_process.nc'
(figure, axes, cbar_axes) = mpl_helper.make_fig(left_margin=0.7, right_margin=1, colorbar=True)

ncf = scipy.io.netcdf_file(filename)
time_grid_edges = ncf.variables["time_grid_edges"].data
entropy_edges = ncf.variables["entropy_edges"].data
time_entropy_dist = ncf.variables["time_entropy_dist"].data

p = axes.pcolor(time_grid_edges, entropy_edges, time_entropy_dist.transpose(),
                norm = matplotlib.colors.LogNorm(), linewidths = 0.1)

axes.set_xscale("linear")
axes.set_xlabel(r"time / h")
axes.set_xlim([0,48])
axes.set_xticks([0, 6, 12, 18, 24, 30, 36, 42, 48])

axes.set_yscale("linear")
axes.set_ylabel(r"particle entropy $H_i$")
axes.set_ylim(0, 3)

axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                       orientation='vertical')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"number conc. $n(t,H_i)$ / $\rm m^{-3}$")

out_filename = "out/urban_plume_time_entropy_all.pdf"
figure.savefig(out_filename)
print out_filename
