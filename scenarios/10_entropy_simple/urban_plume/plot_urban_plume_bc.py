#!/usr/bin/env python

import sys, os
sys.path.append("../../../tool")
import mpl_helper
import scipy.io
import numpy as np
import matplotlib.pyplot
from matplotlib.patches import Polygon
from matplotlib.patches import FancyArrowPatch


(figure, axes_array, cbar_axes_array) \
    = mpl_helper.make_fig_array(1,2,figure_width=8,
                                top_margin=0.45, bottom_margin=0.45,
                                left_margin=1, right_margin=0.65,
                                vert_sep=0.3, horiz_sep=1.5,
                                colorbar="individual",colorbar_location="right")

ncf = scipy.io.netcdf_file("out/urban_plume2_process_bc_new.nc", 'r')
time_grid_edges = ncf.variables["time_grid_edges"].data
entropy_edges = ncf.variables["entropy_edges"].data
time_entropy_dist = ncf.variables["time_entropy_dist"].data


time = ncf.variables["time"].data / 3600
avg_part_entropy = ncf.variables["avg_part_entropy"].data
entropy_of_avg_part = ncf.variables["entropy_of_avg_part"].data
tot_entropy_ratio = ncf.variables["tot_entropy_ratio"].data
ncf.close()

axes = axes_array[0][0]
cbar_axes = cbar_axes_array[0][0]
p = axes.pcolor(time_grid_edges, entropy_edges, time_entropy_dist.transpose(),
                norm = matplotlib.colors.LogNorm(), linewidths = 0.1)

axes.set_xscale("linear")
axes.set_xlabel(r"time / h")
axes.set_xlim([0,48])
axes.set_xticks([0, 6, 12, 18, 24, 30, 36, 42, 48])

axes.set_yscale("linear")
axes.set_ylabel(r"particle entropy $H_i$")
axes.set_ylim(0, 2.5)

axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                       orientation='vertical')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"num. conc. $n(t, H_i)$ / $\rm m^{-3}$")

axes = axes_array[0][1]
axes.plot(time, avg_part_entropy, "b-", label = r'$\bar{H}$')
axes.plot(time, entropy_of_avg_part, "k:", label = r'$\bar{\hat{H}}$')
axes.plot(time, tot_entropy_ratio, "r--", markersize = 2, label = r'$\gamma$')
axes.set_xlabel(r"time / h")
axes.set_ylabel(r"quantity")
axes.set_xlim([0,48])
axes.set_xticks([0, 6, 12, 18, 24, 30, 36, 42, 48])

axes.annotate(r"$\bar{H}$", (time[42], avg_part_entropy[42]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$\bar{\hat{H}}$", (time[66], entropy_of_avg_part[66]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$\gamma$", (time[102], tot_entropy_ratio[102]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')

#axes.legend(loc = 'lower right')
axes.grid(True)

figure.savefig("out/urban_plume_particle_entropy_bc_new.pdf")

