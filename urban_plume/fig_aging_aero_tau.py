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

out_prefix = "figs_aging/aging_aero_tau"

plot_info = {
    "num_low": {"label": r"$\tau_{\rm N}$",
                 "label_time": 12, "label_pos": [0, 0],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[0], "pattern": line_style_list[0],
                 "graph": "g11"},
    "num_low_cond": {"label": r"$\tau^{\rm cond}_{\rm N}$",
                 "label_time": 12, "label_pos": [1, 1],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[2], "pattern": line_style_list[1],
                 "graph": "g11"},
    "num_mid": {"label": r"$\tau_{\rm N}$",
                 "label_time": 12, "label_pos": [0, 0],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[0], "pattern": line_style_list[0],
                 "graph": "g21"},
    "num_mid_cond": {"label": r"$\tau^{\rm cond}_{\rm N}$",
                 "label_time": 12, "label_pos": [1, 1],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[2], "pattern": line_style_list[1],
                 "graph": "g21"},
    "num_high": {"label": r"$\tau_{\rm N}$",
                 "label_time": 12, "label_pos": [0, 0],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[0], "pattern": line_style_list[0],
                 "graph": "g31"},
    "num_high_cond": {"label": r"$\tau^{\rm cond}_{\rm N}$",
                 "label_time": 12, "label_pos": [1, 1],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[2], "pattern": line_style_list[1],
                 "graph": "g31"},
    "mass_low": {"label": r"$\tau_{\rm M}$",
                 "label_time": 12, "label_pos": [0, 0],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[0], "pattern": line_style_list[0],
                 "graph": "g12"},
    "mass_low_cond": {"label": r"$\tau^{\rm cond}_{\rm M}$",
                 "label_time": 12, "label_pos": [1, 1],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[2], "pattern": line_style_list[1],
                 "graph": "g12"},
    "mass_mid": {"label": r"$\tau_{\rm M}$",
                 "label_time": 12, "label_pos": [0, 0],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[0], "pattern": line_style_list[0],
                 "graph": "g22"},
    "mass_mid_cond": {"label": r"$\tau^{\rm cond}_{\rm M}$",
                 "label_time": 12, "label_pos": [1, 1],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[2], "pattern": line_style_list[1],
                 "graph": "g22"},
    "mass_high": {"label": r"$\tau_{\rm M}$",
                 "label_time": 12, "label_pos": [0, 0],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[0], "pattern": line_style_list[0],
                 "graph": "g32"},
    "mass_high_cond": {"label": r"$\tau^{\rm cond}_{\rm M}$",
                 "label_time": 12, "label_pos": [1, 1],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[2], "pattern": line_style_list[1],
                 "graph": "g32"},
    }

time = loadtxt("%s/aging_wc_num_time.txt" % aging_data_dir)
comp_vol = loadtxt("%s/aging_wc_num_comp_vol.txt" % aging_data_dir)

num_tau_transfer = loadtxt("%s/aging_wc_num_tau_transfer.txt" % aging_data_dir)
num_tau_transfer_cond = loadtxt("%s/aging_wc_num_tau_transfer_cond.txt" % aging_data_dir)
mass_tau_transfer = loadtxt("%s/aging_wc_mass_tau_transfer.txt" % aging_data_dir)
mass_tau_transfer_cond = loadtxt("%s/aging_wc_mass_tau_transfer_cond.txt" % aging_data_dir)

num_tau_transfer_smooth = loadtxt("%s/aging_wc_num_tau_transfer_smooth.txt" % aging_data_dir)
num_tau_transfer_cond_smooth = loadtxt("%s/aging_wc_num_tau_transfer_cond_smooth.txt" % aging_data_dir)
mass_tau_transfer_smooth = loadtxt("%s/aging_wc_mass_tau_transfer_smooth.txt" % aging_data_dir)
mass_tau_transfer_cond_smooth = loadtxt("%s/aging_wc_mass_tau_transfer_cond_smooth.txt" % aging_data_dir)

print "low level %d = %f%%" % (level_low, level_low_value * 100)
print "mid level %d = %f%%" % (level_mid, level_mid_value * 100)
print "high level %d = %f%%" % (level_high, level_high_value * 100)

num_low = num_tau_transfer[:,level_low] / 3600 # s to hour
num_mid = num_tau_transfer[:,level_mid] / 3600 # s to hour
num_high = num_tau_transfer[:,level_high] / 3600 # s to hour
num_low_cond = num_tau_transfer_cond[:,level_low] / 3600 # s to hour
num_mid_cond = num_tau_transfer_cond[:,level_mid] / 3600 # s to hour
num_high_cond = num_tau_transfer_cond[:,level_high] / 3600 # s to hour
mass_low = mass_tau_transfer[:,level_low] / 3600 # s to hour
mass_mid = mass_tau_transfer[:,level_mid] / 3600 # s to hour
mass_high = mass_tau_transfer[:,level_high] / 3600 # s to hour
mass_low_cond = mass_tau_transfer_cond[:,level_low] / 3600 # s to hour
mass_mid_cond = mass_tau_transfer_cond[:,level_mid] / 3600 # s to hour
mass_high_cond = mass_tau_transfer_cond[:,level_high] / 3600 # s to hour

num_low_smooth = num_tau_transfer_smooth[:,level_low] / 3600 # s to hour
num_mid_smooth = num_tau_transfer_smooth[:,level_mid] / 3600 # s to hour
num_high_smooth = num_tau_transfer_smooth[:,level_high] / 3600 # s to hour
num_low_cond_smooth = num_tau_transfer_cond_smooth[:,level_low] / 3600 # s to hour
num_mid_cond_smooth = num_tau_transfer_cond_smooth[:,level_mid] / 3600 # s to hour
num_high_cond_smooth = num_tau_transfer_cond_smooth[:,level_high] / 3600 # s to hour
mass_low_smooth = mass_tau_transfer_smooth[:,level_low] / 3600 # s to hour
mass_mid_smooth = mass_tau_transfer_smooth[:,level_mid] / 3600 # s to hour
mass_high_smooth = mass_tau_transfer_smooth[:,level_high] / 3600 # s to hour
mass_low_cond_smooth = mass_tau_transfer_cond_smooth[:,level_low] / 3600 # s to hour
mass_mid_cond_smooth = mass_tau_transfer_cond_smooth[:,level_mid] / 3600 # s to hour
mass_high_cond_smooth = mass_tau_transfer_cond_smooth[:,level_high] / 3600 # s to hour

env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max(time) / 60

for use_color in [True, False]:
    c = canvas.canvas()
    
    x_axis = graph.axis.linear(
        min = 0.,
        max = max_time_min,
        parter = graph.axis.parter.linear(tickdists
                                          = [6 * 60, 3 * 60]),
        texter = time_of_day(base_time
                             = start_time_of_day_min),
        title = "local standard time (LST) (hours:minutes)",
        painter = grid_painter)

    y_axis = graph.axis.log(
        min = 1e-2,
        max = 1e3,
        title = r"aging timescale $\tau$ (hours)",
        painter = major_grid_painter)

    g31 = c.insert(graph.graphxy(
        width = 6.8,
        x = x_axis,
        y = y_axis))
    g21 = c.insert(graph.graphxy(
        width = 6.8,
        ypos = g31.ypos + g31.height + grid_v_space,
        x = graph.axis.linkedaxis(g31.axes["x"],
                                  painter = linked_grid_painter),
        y = y_axis))
    g11 = c.insert(graph.graphxy(
        width = 6.8,
        ypos = g21.ypos + g21.height + grid_v_space,
        x = graph.axis.linkedaxis(g31.axes["x"],
                                  painter = linked_grid_painter),
        y = y_axis))
    g32 = c.insert(graph.graphxy(
        width = 6.8,
        xpos = g31.xpos + g31.width + grid_h_space,
        x = x_axis,
        y = graph.axis.linkedaxis(g31.axes["y"],
                                  painter = linked_major_grid_painter)))
    g22 = c.insert(graph.graphxy(
        width = 6.8,
        xpos = g32.xpos,
        ypos = g21.ypos,
        x = graph.axis.linkedaxis(g32.axes["x"],
                                  painter = linked_grid_painter),
        y = graph.axis.linkedaxis(g21.axes["y"],
                                  painter = linked_major_grid_painter)))
    g12 = c.insert(graph.graphxy(
        width = 6.8,
        xpos = g32.xpos,
        ypos = g11.ypos,
        x = graph.axis.linkedaxis(g32.axes["x"],
                                  painter = linked_grid_painter),
        y = graph.axis.linkedaxis(g11.axes["y"],
                                  painter = linked_major_grid_painter)))

    graphs = {"g11": g11, "g21": g21, "g31": g31,
              "g12": g12, "g22": g22, "g32": g32}

    for (key, y_data) \
            in [("num_low", num_low),
                ("num_low_cond", num_low_cond),
                ("num_mid", num_mid),
                ("num_mid_cond", num_mid_cond),
                ("num_high", num_high),
                ("num_high_cond", num_high_cond),
                ("mass_low", mass_low),
                ("mass_low_cond", mass_low_cond),
                ("mass_mid", mass_mid),
                ("mass_mid_cond", mass_mid_cond),
                ("mass_high", mass_high),
                ("mass_high_cond", mass_high_cond)]:
        g = graphs[plot_info[key]["graph"]]
        if use_color:
            grey_color = color.hsb(plot_info[key]["color"].hsb().color["h"], grey_level, 1)
            style_attrs = [plot_info[key]["linewidth"],
                           grey_color]
        else:
            grey_color = color.grey(1 - grey_level)
            style_attrs = [plot_info[key]["linewidth"],
                           grey_color]
        plot_data = zip(time[1:] / 60, y_data)
        g.plot(
            graph.data.points(plot_data, x = 1, y = 2),
            styles = [graph.style.line(lineattrs = style_attrs)])

    for (key, y_data) \
            in [("num_low", num_low_smooth),
                ("num_low_cond", num_low_cond_smooth),
                ("num_mid", num_mid_smooth),
                ("num_mid_cond", num_mid_cond_smooth),
                ("num_high", num_high_smooth),
                ("num_high_cond", num_high_cond_smooth),
                ("mass_low", mass_low_smooth),
                ("mass_low_cond", mass_low_cond_smooth),
                ("mass_mid", mass_mid_smooth),
                ("mass_mid_cond", mass_mid_cond_smooth),
                ("mass_high", mass_high_smooth),
                ("mass_high_cond", mass_high_cond_smooth)]:
        g = graphs[plot_info[key]["graph"]]
        if use_color:
            style_attrs = [plot_info[key]["linewidth"],
                           plot_info[key]["color"]]
        else:
            style_attrs = [plot_info[key]["linewidth"],
                           plot_info[key]["pattern"]]
        plot_data = zip(time[1:] / 60, y_data)
        g.plot(
            graph.data.points(plot_data, x = 1, y = 2),
            styles = [graph.style.line(lineattrs = style_attrs)])

    for (g_name, g) in graphs.iteritems():
        g.dodata()
        g.doaxes()

    for (key, y_data) \
            in [("num_low", num_low_smooth),
                ("num_low_cond", num_low_cond_smooth),
                ("num_mid", num_mid_smooth),
                ("num_mid_cond", num_mid_cond_smooth),
                ("num_high", num_high_smooth),
                ("num_high_cond", num_high_cond_smooth),
                ("mass_low", mass_low_smooth),
                ("mass_low_cond", mass_low_cond_smooth),
                ("mass_mid", mass_mid_smooth),
                ("mass_mid_cond", mass_mid_cond_smooth),
                ("mass_high", mass_high_smooth),
                ("mass_high_cond", mass_high_cond_smooth)]:
        g = graphs[plot_info[key]["graph"]]
        plot_data = zip(time[1:] / 60, y_data)
        label_plot_line_boxed(g, plot_data,
                              plot_info[key]["label_time"] * 60,
                              plot_info[key]["label"],
                              plot_info[key]["label_pos"])

    write_text_outside(g11, r"critical supersaturation $S = %.1f\%%$" % (level_low_value * 100))
    write_text_outside(g21, r"critical supersaturation $S = %.1f\%%$" % (level_mid_value * 100))
    write_text_outside(g31, r"critical supersaturation $S = %.1f\%%$" % (level_high_value * 100))
    write_text_outside(g12, r"critical supersaturation $S = %.1f\%%$" % (level_low_value * 100))
    write_text_outside(g22, r"critical supersaturation $S = %.1f\%%$" % (level_mid_value * 100))
    write_text_outside(g32, r"critical supersaturation $S = %.1f\%%$" % (level_high_value * 100))

    boxed_text(g11, "number", point = [1, 1], anchor_point_rel = [1, 1])
    boxed_text(g21, "number", point = [1, 1], anchor_point_rel = [1, 1])
    boxed_text(g31, "number", point = [1, 1], anchor_point_rel = [1, 1])
    boxed_text(g12, "mass", point = [1, 1], anchor_point_rel = [1, 1])
    boxed_text(g22, "mass", point = [1, 1], anchor_point_rel = [1, 1])
    boxed_text(g32, "mass", point = [1, 1], anchor_point_rel = [1, 1])

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    c.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
