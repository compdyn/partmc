#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import matplotlib
import partmc
import scipy.io, numpy

(figure, axes_array, cbar_axes_array) \
    = mpl_helper.make_fig_array(4,2,figure_width=8,
                                top_margin=0.45, bottom_margin=0.45,
                                left_margin=1, right_margin=0.65,
                                vert_sep=0.3, horiz_sep=1.5,
                                colorbar="individual",colorbar_location="right",
                                share_y_axes=False)

######## first row ###############
filename = 'coag/out/urban_plume_process.nc'
ncf = scipy.io.netcdf_file(filename)
time_grid_edges = ncf.variables["time_grid_edges"].data
diversity_edges = ncf.variables["diversity_edges"].data
time_diversity_dist = ncf.variables["time_diversity_dist"].data

time = ncf.variables["time"].data / 3600
avg_part_entropy = numpy.exp(ncf.variables["avg_part_entropy"].data)
entropy_of_avg_part = numpy.exp(ncf.variables["entropy_of_avg_part"].data)
#tot_entropy_ratio = ncf.variables["tot_entropy_ratio"].data
tot_entropy_ratio = (avg_part_entropy-1) / (entropy_of_avg_part-1)
ncf.close()

axes = axes_array[0][0]
cbar_axes = cbar_axes_array[0][0]
p = axes.pcolor(time_grid_edges, diversity_edges, time_diversity_dist.transpose(),
                norm = matplotlib.colors.LogNorm(), linewidths = 0.1)

axes.set_xscale("linear")
axes.set_xlabel(r"time / h")
axes.set_xlim([0,24])
axes.set_xticks([0, 6, 12, 18, 24])

axes.set_yscale("linear")
axes.set_ylabel(r"particle diversity $D_i$")
axes.set_ylim(0.95, 2.05)

axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                       orientation='vertical')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"num. conc. $n(t, D_i)$ / $\rm m^{-3}$")

axes = axes_array[0][1]
axes.plot(time, avg_part_entropy, "b-")
axes.plot(time, entropy_of_avg_part, "k:")
axes.set_xlabel(r"time / h")
axes.set_xlim([0,24])
axes.set_xticks([0, 6, 12, 18, 24])

axes.set_ylabel(r"$D_{\alpha}, D{\gamma}$")
axes.set_ylim([0.95,2.05])

axes2 =  axes.twinx()
axes2.plot(time, tot_entropy_ratio, "r--", markersize = 2)
axes2.set_ylabel(r"$\chi$")

axes.annotate(r"$D_{\alpha}$", (time[78], avg_part_entropy[78]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$D_{\gamma}$", (time[66], entropy_of_avg_part[66]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$\chi$", (time[102], tot_entropy_ratio[102]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')

axes.grid(True)

######## second row ###############
filename = 'emit_pure/out/urban_plume_process.nc'
ncf = scipy.io.netcdf_file(filename)
time_grid_edges = ncf.variables["time_grid_edges"].data
diversity_edges = ncf.variables["diversity_edges"].data
time_diversity_dist = ncf.variables["time_diversity_dist"].data

time = ncf.variables["time"].data / 3600
avg_part_entropy = numpy.exp(ncf.variables["avg_part_entropy"].data)
entropy_of_avg_part = numpy.exp(ncf.variables["entropy_of_avg_part"].data)
#tot_entropy_ratio = ncf.variables["tot_entropy_ratio"].data
tot_entropy_ratio = (avg_part_entropy-1) / (entropy_of_avg_part-1)
ncf.close()

axes = axes_array[1][0]
cbar_axes = cbar_axes_array[1][0]
p = axes.pcolor(time_grid_edges, diversity_edges, time_diversity_dist.transpose(),
                norm = matplotlib.colors.LogNorm(), linewidths = 0.1)

axes.set_xscale("linear")
axes.set_xlabel(r"time / h")
axes.set_xlim([0,24])
axes.set_xticks([0, 6, 12, 18, 24])

axes.set_yscale("linear")
axes.set_ylabel(r"particle diversity $D_i$")
axes.set_ylim(0.95, 2.05)

axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                       orientation='vertical')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"num. conc. $n(t, D_i)$ / $\rm m^{-3}$")

axes = axes_array[1][1]
axes.plot(time, avg_part_entropy, "b-")
axes.plot(time, entropy_of_avg_part, "k:")
axes.set_xlabel(r"time / h")
axes.set_xlim([0,24])
axes.set_xticks([0, 6, 12, 18, 24])

axes.set_ylabel(r"$D_{\alpha}, D{\gamma}$")
axes.set_ylim([0.95,3.0])

axes2 =  axes.twinx()
axes2.plot(time, tot_entropy_ratio, "r--", markersize = 2)
axes2.set_ylabel(r"$\chi$")

axes.annotate(r"$D_{\alpha}$", (time[78], avg_part_entropy[78]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$D_{\gamma}$", (time[66], entropy_of_avg_part[66]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$\chi$", (time[102], tot_entropy_ratio[102]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')

axes.grid(True)

######## third row ###############
filename = 'emit_mixed/out/urban_plume_process.nc'
ncf = scipy.io.netcdf_file(filename)
time_grid_edges = ncf.variables["time_grid_edges"].data
diversity_edges = ncf.variables["diversity_edges"].data
time_diversity_dist = ncf.variables["time_diversity_dist"].data

time = ncf.variables["time"].data / 3600
avg_part_entropy = numpy.exp(ncf.variables["avg_part_entropy"].data)
entropy_of_avg_part = numpy.exp(ncf.variables["entropy_of_avg_part"].data)
#tot_entropy_ratio = ncf.variables["tot_entropy_ratio"].data
tot_entropy_ratio = (avg_part_entropy-1) / (entropy_of_avg_part-1)
ncf.close()

axes = axes_array[2][0]
cbar_axes = cbar_axes_array[2][0]
p = axes.pcolor(time_grid_edges, diversity_edges, time_diversity_dist.transpose(),
                norm = matplotlib.colors.LogNorm(), linewidths = 0.1)

axes.set_xscale("linear")
axes.set_xlabel(r"time / h")
axes.set_xlim([0,24])
axes.set_xticks([0, 6, 12, 18, 24])

axes.set_yscale("linear")
axes.set_ylabel(r"particle diversity $D_i$")
axes.set_ylim(1, 2.5)

axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                       orientation='vertical')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"num. conc. $n(t, D_i)$ / $\rm m^{-3}$")

axes = axes_array[2][1]
axes.plot(time, avg_part_entropy, "b-")
axes.plot(time, entropy_of_avg_part, "k:")
axes.set_xlabel(r"time / h")
axes.set_xlim([0,24])
axes.set_xticks([0, 6, 12, 18, 24])

axes.set_ylabel(r"$D_{\alpha}, D{\gamma}$")
axes.set_ylim([1,5.0])

axes2 =  axes.twinx()
axes2.plot(time, tot_entropy_ratio, "r--", markersize = 2)
axes2.set_ylabel(r"$\chi$")

axes.annotate(r"$D_{\alpha}$", (time[78], avg_part_entropy[78]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$D_{\gamma}$", (time[66], entropy_of_avg_part[66]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$\chi$", (time[102], tot_entropy_ratio[102]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')

axes.grid(True)

######## forth row ###############
filename = 'emit/out/urban_plume_process.nc'
ncf = scipy.io.netcdf_file(filename)
time_grid_edges = ncf.variables["time_grid_edges"].data
diversity_edges = ncf.variables["diversity_edges"].data
time_diversity_dist = ncf.variables["time_diversity_dist"].data

time = ncf.variables["time"].data / 3600
avg_part_entropy = numpy.exp(ncf.variables["avg_part_entropy"].data)
entropy_of_avg_part = numpy.exp(ncf.variables["entropy_of_avg_part"].data)
#tot_entropy_ratio = ncf.variables["tot_entropy_ratio"].data
tot_entropy_ratio = (avg_part_entropy-1) / (entropy_of_avg_part-1)
ncf.close()

axes = axes_array[3][0]
cbar_axes = cbar_axes_array[3][0]
p = axes.pcolor(time_grid_edges, diversity_edges, time_diversity_dist.transpose(),
                norm = matplotlib.colors.LogNorm(), linewidths = 0.1)

axes.set_xscale("linear")
axes.set_xlabel(r"time / h")
axes.set_xlim([0,24])
axes.set_xticks([0, 6, 12, 18, 24])

axes.set_yscale("linear")
axes.set_ylabel(r"particle diversity $D_i$")
axes.set_ylim(1, 2.5)

axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                       orientation='vertical')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"num. conc. $n(t, D_i)$ / $\rm m^{-3}$")

axes = axes_array[3][1]
axes.plot(time, avg_part_entropy, "b-")
axes.plot(time, entropy_of_avg_part, "k:")
axes.set_xlabel(r"time / h")
axes.set_xlim([0,24])
axes.set_xticks([0, 6, 12, 18, 24])

axes.set_ylabel(r"$D_{\alpha}, D{\gamma}$")
axes.set_ylim([1,5.0])

axes2 =  axes.twinx()
axes2.plot(time, tot_entropy_ratio, "r--", markersize = 2)
axes2.set_ylabel(r"$\chi$")

axes.annotate(r"$D_{\alpha}$", (time[78], avg_part_entropy[78]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$D_{\gamma}$", (time[66], entropy_of_avg_part[66]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$\chi$", (time[102], tot_entropy_ratio[102]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')

axes.grid(True)

#mpl_helper.remove_fig_array_axes(axes_array)

out_filename = "simple_cases1.pdf"
figure.savefig(out_filename)
print out_filename
