#!/usr/bin/env python

import sys, os
sys.path.append("../../../tool")
import mpl_helper
import scipy.io, numpy

(figure, axes) = mpl_helper.make_fig(right_margin=0.8)

ncf = scipy.io.netcdf_file("out/urban_plume_process.nc")
time = ncf.variables["time"].data / 3600

avg_part_entropy = ncf.variables["avg_part_entropy"].data
avg_part_entropy_err = ncf.variables["avg_part_entropy_ci_offset"].data
entropy_of_avg_part = ncf.variables["entropy_of_avg_part"].data
entropy_of_avg_part_err = ncf.variables["entropy_of_avg_part_ci_offset"].data
tot_entropy_ratio = ncf.variables["tot_entropy_ratio"].data

axes.plot(time, avg_part_entropy, "b-", label = r'$\bar{H}$')
axes.plot(time, entropy_of_avg_part, "k:", label = r'$\bar{\hat{H}}$')
axes.plot(time, tot_entropy_ratio, "r--", markersize = 2, label = r'$\gamma$')
axes.set_xlabel(r"time / h")
axes.set_ylabel(r"quantity")
axes.set_xlim([0,24])
axes.set_xticks([0, 6, 12, 18, 24])

axes.annotate(r"$\bar{H}$", (time[7], avg_part_entropy[7]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$\bar{\hat{H}}$", (time[11], entropy_of_avg_part[11]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$\gamma$", (time[17], tot_entropy_ratio[17]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')

#axes.legend(loc = 'lower right')
axes.grid(True)

figure.savefig("out/urban_plume_particle_entropy.pdf")
