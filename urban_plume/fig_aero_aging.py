#!/usr/bin/env python
# Copyright (C) 2007, 2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
from fig_helper import *

out_prefix = "figs/aero_aging"

time_filename_list_wc = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max([time for [time, filename] in time_filename_list_wc]) / 60

g = graph.graphxy(
    width = 6.7,
    x = graph.axis.linear(min = 0.,
                          max = max_time_min,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time
                                               = start_time_of_day_min),
                          title = "local standard time (hours:minutes)",
                          painter = grid_painter),
    y = graph.axis.linear(title = r"mass density ($\rm \mu g \, m^{-3}$)",
                          painter = grid_painter),
    y2 = graph.axis.linear(title = r"number density ($\rm m^{-3}$)"))

plot_data = [[] for i in range(len(plots))]
for [time, filename] in time_filename_list_wc:
    print i, filename
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    ncf.close()
    value = (particles.mass() / particles.comp_vol).sum() * 1e9
    value = (particles.mass(exclude = ["H2O"])
             / particles.comp_vol).sum() * 1e9
    value = (1.0 / particles.comp_vol).sum()
    plot_data[i].append([time / 60.0, value])

for i in range(len(plots)):
    attrs = []
    if use_color:
        attrs.append(color_list[plots[i]["color_index"]])
        attrs.append(plots[i]["color_style"])
    else:
        attrs.append(line_style_list[i])
    attrs.append(style.linewidth.Thick)
    if plots[i]["mass"]:
        g.plot(
            graph.data.points(
            plot_data[i], x = 1, y = 2),
            styles = [graph.style.line(lineattrs = attrs)])
    else:
        g.plot(
            graph.data.points(
            plot_data[i], x = 1, y2 = 2),
            styles = [graph.style.line(lineattrs = attrs)])

for i in range(len(plots)):
    if plots[i]["label"] != None:
        label_plot_line(g, plot_data[i], plots[i]["time"] * 60.0,
                        plots[i]["label"], plots[i]["pos"],
                        yaxis = g.axes[plots[i]["yaxis"]])

if not use_color:
    for i in range(len(plots)):
        if (plots[i]["mass"] == False):
            peak_val = max([v for [t, v] in plot_data[i]])
            final_val = plot_data[i][-1][1]
            print "total number, coag = %s, peak = %g, final = %g m^{-3}" \
                  % (str(plots[i]["coag"]), peak_val, final_val)
            if plots[i]["coag"] == False:
                base_peak = peak_val
                base_final = final_val
            else:
                new_peak = peak_val
                new_final = final_val
    print "total number, peak_val decrease from no-coag to with-coag = %g%%" \
          % ((base_peak - new_peak) / base_peak * 100)
    print "total number, final_val decrease from no-coag to with-coag = %g%%" \
          % ((base_final - new_final) / base_final * 100)

if use_color:
    out_filename = "%s_color.pdf" % out_prefix
else:
    out_filename = "%s_bw.pdf" % out_prefix
g.writePDFfile(out_filename)
if not use_color:
    print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
    print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
