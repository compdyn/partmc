#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import matplotlib
import partmc
import scipy.io, numpy as np

(figure, axes_array) \
    = mpl_helper.make_fig_array(3,1,figure_width=4,
                                top_margin=0.1, bottom_margin=0.45,
                                left_margin=0.5, right_margin=1,
                                vert_sep=0.2,
                                share_y_axes=False)

######## first row ###############
filename = 'coag/out/urban_plume_process.nc'
ncf = scipy.io.netcdf_file(filename)

time = ncf.variables["time"].data / 3600
avg_part_entropy = np.exp(ncf.variables["avg_part_entropy"].data)
entropy_of_avg_part = np.exp(ncf.variables["entropy_of_avg_part"].data)
#tot_entropy_ratio = ncf.variables["tot_entropy_ratio"].data
tot_entropy_ratio = (avg_part_entropy-1) / (entropy_of_avg_part-1)
ncf.close()

axes = axes_array[0][0]
axes.plot(time, avg_part_entropy, "b-")
axes.plot(time, entropy_of_avg_part, "k--")

axes.set_ylabel(r"$D_{\alpha}, D{\gamma}$")
axes.set_ylim([0.95,2.05])

axes2 =  axes.twinx()
axes2.plot(time, tot_entropy_ratio, "r-", markersize = 2)
axes2.set_ylabel(r"$\chi$")
axes2.set_yticks([0,0.2,0.4,0.6,0.8])

axes2.set_xscale("linear")
axes2.set_xlabel(r"time / h")
axes2.set_xlim([0,24])
axes2.set_xticks([0, 6, 12, 18, 24])

print "coordinates "
print time[78], avg_part_entropy[78]
print time[66], entropy_of_avg_part[66]
print time[102], tot_entropy_ratio[102]

axes.annotate(r"$D_{\alpha}$", (20,1.5),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$D_{\gamma}$", (3,1.8),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$\chi$", (10,1.7),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')

axes.annotate(r"Coagulation", (2,1.05),
              verticalalignment="bottom", horizontalalignment="left",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')

axes.grid(True)

######## second row ###############
filename = 'emit/out/urban_plume_process.nc'
ncf = scipy.io.netcdf_file(filename)

time = ncf.variables["time"].data / 3600
avg_part_entropy = np.exp(ncf.variables["avg_part_entropy"].data)
entropy_of_avg_part = np.exp(ncf.variables["entropy_of_avg_part"].data)
#tot_entropy_ratio = ncf.variables["tot_entropy_ratio"].data
tot_entropy_ratio = (avg_part_entropy-1) / (entropy_of_avg_part-1)
ncf.close()

axes = axes_array[1][0]
axes.plot(time, avg_part_entropy, "b-")
axes.plot(time, entropy_of_avg_part, "k--")
axes.set_ylabel(r"$D_{\alpha}, D{\gamma}$")
axes.set_ylim([0.95,3.0])

axes2 =  axes.twinx()
axes2.plot(time, tot_entropy_ratio, "r-", markersize = 2)
axes2.set_ylabel(r"$\chi$")
axes2.set_xscale("linear")
axes2.set_xlim([0,24])
axes2.set_xticks([0, 6, 12, 18, 24])


axes.annotate(r"$D_{\alpha}$", (time[78], avg_part_entropy[78]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$D_{\gamma}$", (time[66], entropy_of_avg_part[66]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$\chi$", (20,1.15),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')

axes.annotate(r"Emission", (2,1.1),
              verticalalignment="bottom", horizontalalignment="left",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.grid(True)

######## third row ###############
filename = 'cond_mono/out/urban_plume_process.nc'
ncf = scipy.io.netcdf_file(filename)

time = ncf.variables["time"].data / 3600
avg_part_entropy = np.exp(ncf.variables["avg_part_entropy"].data)
entropy_of_avg_part = np.exp(ncf.variables["entropy_of_avg_part"].data)
#tot_entropy_ratio = ncf.variables["tot_entropy_ratio"].data
tot_entropy_ratio = (avg_part_entropy-1) / (entropy_of_avg_part-1)
ncf.close()

axes = axes_array[2][0]
axes.plot(time, avg_part_entropy, "b-")
axes.plot(time, entropy_of_avg_part, "k--")

axes.set_ylabel(r"$D_{\alpha}, D{\gamma}$")
axes.set_ylim([0.95,2.5])
axes.set_yticks([1, 1.5, 2, 2.5])

axes2 =  axes.twinx()
axes2.plot(time, tot_entropy_ratio, "r-", markersize = 2)
axes2.set_ylabel(r"$\chi$")
axes2.set_xscale("linear")
axes2.set_xlim([0,24])
axes2.set_xticks([0, 6, 12, 18, 24])

axes.annotate(r"$D_{\alpha}$", (10,2.27),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$D_{\gamma}$", (8,1.95),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$\chi$", (18,1.7),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"Condensation", (2,1.1),
              verticalalignment="bottom", horizontalalignment="left",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')

axes.grid(True)

mpl_helper.remove_fig_array_axes(axes_array)

out_filename = "simple_poster2.pdf"
figure.savefig(out_filename)
print out_filename
