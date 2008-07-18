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

out_filename = "figs/aero_totals.pdf"

plots = [
    {"mass": True, "coag": True, "wet": True, "yaxis": "y",
     "label": "wet mass", "time": 10, "pos": [1, 1]},
    {"mass": True, "coag": False, "wet": True, "yaxis": "y",
     "label": None, "time": 10, "pos": [0, 0]},
    {"mass": False, "coag": True, "yaxis": "y2",
     "label": "num with coag", "time": 12, "pos": [0.3, 1]},
    {"mass": False, "coag": False, "yaxis": "y2",
     "label": "num no coag", "time": 8, "pos": [1, 1]},
    {"mass": True, "coag": True, "wet": False, "yaxis": "y",
     "label": "dry mass", "time": 15, "pos": [1, 0]},
    {"mass": True, "coag": False, "wet": False, "yaxis": "y",
     "label": None, "time": 18, "pos": [0, 0]},
    ]

time_filename_list_wc = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
time_filename_list_nc = get_time_filename_list(netcdf_dir_nc, netcdf_pattern_nc)
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
for i in range(len(plots)):
    if plots[i]["coag"]:
        time_filename_list = time_filename_list_wc
    else:
        time_filename_list = time_filename_list_nc
    for [time, filename] in time_filename_list:
        print i, filename
        ncf = NetCDFFile(filename)
        particles = aero_particle_array_t(ncf)
        ncf.close()
        if plots[i]["mass"]:
            if plots[i]["wet"]:
                value = (particles.mass() / particles.comp_vol).sum() * 1e9
            else:
                value = (particles.mass(exclude = ["H2O"])
                         / particles.comp_vol).sum() * 1e9
        else:
            value = (1.0 / particles.comp_vol).sum()
        plot_data[i].append([time / 60.0, value])

for i in range(len(plots)):
    if plots[i]["mass"]:
        g.plot(
            graph.data.points(
            plot_data[i], x = 1, y = 2),
            styles = [graph.style.line(lineattrs = [line_style_list[i],
                                                    style.linewidth.Thick])])
    else:
        g.plot(
            graph.data.points(
            plot_data[i], x = 1, y2 = 2),
            styles = [graph.style.line(lineattrs = [line_style_list[i],
                                                    style.linewidth.Thick])])

for i in range(len(plots)):
    if plots[i]["label"] != None:
        label_plot_line(g, plot_data[i], plots[i]["time"] * 60.0,
                        plots[i]["label"], plots[i]["pos"],
                        yaxis = g.axes[plots[i]["yaxis"]])

for i in range(len(plots)):
    if (plots[i]["mass"] == False):
        print "total number, coag = %s, max = %g, final = %g m^{-3}" \
              % (str(plots[i]["coag"]), max([v for [t, v] in plot_data[i]]),
                 plot_data[i][-1][1])

g.writePDFfile(out_filename)
print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
