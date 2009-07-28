#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
sys.path.append(".")
from fig_helper import *

data_file = "out/aging_aero_time_bulk.txt"
out_prefix = "figs_aging/aging_aero_time_bulk"

plot_info = {
    "num": {
        "label": r"number", "label_time": 8,
        "label_xoffset": 0 * unit.v_mm, "label_yoffset": 0 * unit.v_mm,
        "label_pos": [1, 1],
        "linewidth": style.linewidth.Thick,
        "color": color_list[0], "pattern": line_style_list[0]},
    "mass": {
        "label": r"wet mass", "label_time": 11.5,
        "label_xoffset": 0 * unit.v_mm, "label_yoffset": 0 * unit.v_mm,
        "label_pos": [1, 1],
        "linewidth": style.linewidth.Thick,
        "color": color_list[1], "pattern": line_style_list[1]},
    "dry_mass": {
        "label": r"dry mass", "label_time": 18.5,
        "label_xoffset": 0 * unit.v_mm, "label_yoffset": 0 * unit.v_mm,
        "label_pos": [0, 1],
        "linewidth": style.linewidth.Thick,
        "color": color_list[1], "pattern": line_style_list[2]},
    "area": {
        "label": r"wet area", "label_time": 18.5,
        "label_xoffset": 0 * unit.v_mm, "label_yoffset": -1 * unit.v_mm,
        "label_pos": [0, 1],
        "linewidth": style.linewidth.Thick,
        "color": color_list[2], "pattern": line_style_list[3]},
    "dry_area": {
        "label": r"dry area", "label_time": 17.5,
        "label_xoffset": 0 * unit.v_mm, "label_yoffset": 0 * unit.v_mm,
        "label_pos": [1, 0],
        "linewidth": style.linewidth.Thick,
        "color": color_list[2], "pattern": line_style_list[4]},
    }

plot_data = {}
for quantity in plot_info.keys():
    plot_data[quantity] = []
data = loadtxt(data_file)
for i_time in range(size(data,0)):
    time_min = data[i_time, 0] / 60.0
    plot_data["num"].append([time_min, data[i_time, 1] * 1e-6])
    plot_data["mass"].append([time_min, data[i_time, 2] * 1e9])
    plot_data["dry_mass"].append([time_min, data[i_time, 3] * 1e9])
    plot_data["area"].append([time_min, data[i_time, 4] * 1e6])
    plot_data["dry_area"].append([time_min, data[i_time, 5] * 1e6])

env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max(data[:,0]) / 60.0

num_values = [n for [t,n] in plot_data["num"]]
print "max number density = %f cm^{-3}" % max(num_values)

for use_color in [True, False]:
    g = graph.graphxy(
        width = 6.8,
        x = graph.axis.linear(min = 0.,
                              max = max_time_min,
                              parter = graph.axis.parter.linear(tickdists
                                                                = [6 * 60, 3 * 60]),
                              texter = time_of_day(base_time
                                                   = start_time_of_day_min),
                              title = "local standard time (LST) (hours:minutes)",
                              painter = grid_painter),
        y = graph.axis.linear(min = 0,
                              max = 2e4,
                              title = r"number conc. ($\rm cm^{-3}$)",
                              painter = grid_painter),
        y2 = graph.axis.linear(min = 0,
                              max = 200,
                              title = r"mass conc. ($\rm \mu g \, m^{-3}$)",
                              painter = grid_painter),
        y4 = graph.axis.linear(min = 0,
                              max = 4e3,
                              title = r"surface area conc. ($\rm \mu m^2 \, cm^{-3}$)",
                              painter = grid_painter))

    g.doaxes()

    y_axes = {"num": "y",
              "mass": "y2",
              "dry_mass": "y2",
              "area": "y4",
              "dry_area": "y4"}
    for quantity in plot_info.keys():
        if use_color:
            style_attrs = [plot_info[quantity]["linewidth"],
                           plot_info[quantity]["color"]]
        else:
            style_attrs = [plot_info[quantity]["linewidth"],
                           plot_info[quantity]["pattern"]]
        y_arg = {y_axes[quantity]: 2}
        g.plot(
            graph.data.points(plot_data[quantity], x = 1, **y_arg),
            styles = [graph.style.line(lineattrs = style_attrs)])
        label_plot_line_boxed(g, plot_data[quantity],
                              plot_info[quantity]["label_time"] * 60,
                              plot_info[quantity]["label"],
                              plot_info[quantity]["label_pos"],
                              label_xoffset = plot_info[quantity]["label_xoffset"],
                              label_yoffset = plot_info[quantity]["label_yoffset"],
                              yaxis = g.axes[y_axes[quantity]])

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    g.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
