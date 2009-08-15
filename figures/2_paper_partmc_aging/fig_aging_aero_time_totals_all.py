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

data_file = "out/aging_aero_time_totals_all.txt"
out_prefix = "figs_aging/aging_aero_time_totals_all"

plot_info = {
    "num" : {
        "total": {"label": r"$N_{\rm BC}$", "label_time": 8,
                  "label_pos": [1, 1],
                  "linewidth": style.linewidth.Thick,
                  "color": color_list[0], "pattern": line_style_list[0]},
        "aged":  {"label": r"$N_{\rm a}$", "label_time": 11.5,
                  "label_pos": [0, 1],
                  "linewidth": style.linewidth.Thick,
                  "color": color_list[1], "pattern": line_style_list[1]},
        "fresh": {"label": r"$N_{\rm f}$", "label_time": 13.5,
                  "label_pos": [1, 0],
                  "linewidth": style.linewidth.Thick,
                  "color": color_list[2], "pattern": line_style_list[2]},
        },
    "mass" : {
        "total": {"label": r"$M^{\rm total}_{\rm BC}$", "label_time": 8,
                  "label_pos": [1, 1],
                  "linewidth": style.linewidth.Thick,
                  "color": color_list[0], "pattern": line_style_list[0]},
        "aged":  {"label": r"$M^{\rm total}_{\rm a}$", "label_time": 11.5,
                  "label_pos": [0, 0],
                  "linewidth": style.linewidth.Thick,
                  "color": color_list[1], "pattern": line_style_list[1]},
        "fresh": {"label": r"$M^{\rm total}_{\rm f}$", "label_time": 13.5,
                  "label_pos": [0, 1],
                  "linewidth": style.linewidth.Thick,
                  "color": color_list[2], "pattern": line_style_list[2]},
        },
    "area" : {
        "total": {"label": r"$A_{\rm BC}$", "label_time": 8,
                  "label_pos": [1, 1],
                  "linewidth": style.linewidth.Thick,
                  "color": color_list[0], "pattern": line_style_list[0]},
        "aged":  {"label": r"$A_{\rm a}$", "label_time": 11.5,
                  "label_pos": [0, 0],
                  "linewidth": style.linewidth.Thick,
                  "color": color_list[1], "pattern": line_style_list[1]},
        "fresh": {"label": r"$A_{\rm f}$", "label_time": 13.5,
                  "label_pos": [0, 1],
                  "linewidth": style.linewidth.Thick,
                  "color": color_list[2], "pattern": line_style_list[2]},
        },
    }

plot_data = {}
for quantity in plot_info.keys():
    plot_data[quantity] = {}
    for age in plot_info[quantity].keys():
        plot_data[quantity][age] = []
data = loadtxt(data_file)
for i_time in range(size(data,0)):
    time_min = data[i_time, 0] / 60.0
    plot_data["num"]["total"].append([time_min, data[i_time, 1] * 1e-6])
    plot_data["num"]["aged"].append([time_min, data[i_time, 2] * 1e-6])
    plot_data["num"]["fresh"].append([time_min, data[i_time, 3] * 1e-6])

    plot_data["mass"]["total"].append([time_min, data[i_time, 4] * 1e9])
    plot_data["mass"]["aged"].append([time_min, data[i_time, 5] * 1e9])
    plot_data["mass"]["fresh"].append([time_min, data[i_time, 6] * 1e9])

    plot_data["area"]["total"].append([time_min, data[i_time, 7] * 1e6])
    plot_data["area"]["aged"].append([time_min, data[i_time, 8] * 1e6])
    plot_data["area"]["fresh"].append([time_min, data[i_time, 9] * 1e6])

env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max(data[:,0]) / 60.0

for use_color in [True, False]:
    c = canvas.canvas()

    g3 = c.insert(graph.graphxy(
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
                              max = 2000,
                              #density = 1.5, # hack to increase number of ticks
                              title = r"surface area conc. ($\rm \mu m^2 cm^{-3}$)",
                              painter = grid_painter)))
    g2 = c.insert(graph.graphxy(
        width = g3.width,
        ypos = g3.ypos + g3.height + 0.5,
        x = graph.axis.linkedaxis(g3.axes["x"],
                                  painter = graph.axis.painter.linked(gridattrs = [style.linestyle.dotted])),
        y = graph.axis.linear(min = 0.0,
                              max = 200.0,
                              title = r"mass conc. ($\rm \mu g \, m^{-3}$)",
                              painter = grid_painter)))
    g1 = c.insert(graph.graphxy(
        width = g2.width,
        ypos = g2.ypos + g2.height + 0.5,
        x = graph.axis.linkedaxis(g3.axes["x"],
                                  painter = graph.axis.painter.linked(gridattrs = [style.linestyle.dotted])),
        y = graph.axis.linear(min = 0.0,
                              max = 1.5e4,
                              title = r"number conc. ($\rm cm^{-3}$)",
                              painter = grid_painter)))

    g1.doaxes()
    g2.doaxes()
    g3.doaxes()

    graphs_by_quantity = {"num": g1,
                          "mass": g2,
                          "area": g3}
    for quantity in plot_info.keys():
        g = graphs_by_quantity[quantity]
        for age in plot_info[quantity].keys():
            if use_color:
                style_attrs = [plot_info[quantity][age]["linewidth"],
                               plot_info[quantity][age]["color"]]
            else:
                style_attrs = [plot_info[quantity][age]["linewidth"],
                               plot_info[quantity][age]["pattern"]]
            g.plot(
                graph.data.points(plot_data[quantity][age], x = 1, y = 2),
                styles = [graph.style.line(lineattrs = style_attrs)])
            label_plot_line_boxed(g, plot_data[quantity][age],
                                  plot_info[quantity][age]["label_time"] * 60,
                                  plot_info[quantity][age]["label"],
                                  plot_info[quantity][age]["label_pos"])

    write_text_outside(g1, r"critical supersaturation $S_{\rm c} = %.1f\%%$" % (level_mid_value * 100))

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    c.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
