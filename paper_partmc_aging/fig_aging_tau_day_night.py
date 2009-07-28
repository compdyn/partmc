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

out_prefix = "figs_aging/aging_tau_day_night"

plot_info = {
    "day_num": {"label": r"$\tau_{\rm N,day}$",
                 "label_time": 0.32, "label_pos": [0, 1],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[0], "pattern": line_style_list[0],
                 "graph": "g1"},
#    "day_num_cond": {"label": r"$\tau^{\rm cond}_{\rm N,day}$",
#                 "label_time": 0.6, "label_pos": [0, 1],
#                 "linewidth": style.linewidth.Thick,
#                 "color": color_list[1], "pattern": line_style_list[1],
#                 "graph": "g1"},
    "night_num": {"label": r"$\tau_{\rm N,night}$",
                 "label_time": 0.065, "label_pos": [1, 0],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[0], "pattern": line_style_list[0],
                 "graph": "g2"},
    "night_num_cond": {"label": r"$\tau^{\rm cond}_{\rm N,night}$",
                 "label_time": 0.65, "label_pos": [0, 1],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[1], "pattern": line_style_list[2],
                 "graph": "g2"},
    "day_mass": {"label": r"$\tau_{\rm M,day}$",
                 "label_time": 0.35, "label_pos": [1, 0],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[2], "pattern": line_style_list[1],
                 "graph": "g1"},
#    "day_mass_cond": {"label": r"$\tau^{\rm cond}_{\rm M,day}$",
#                 "label_time": 0.6, "label_pos": [0, 1],
#                 "linewidth": style.linewidth.Thick,
#                 "color": color_list[3], "pattern": line_style_list[3],
#                 "graph": "g1"},
    "night_mass": {"label": r"$\tau_{\rm M,night}$",
                 "label_time": 0.065, "label_pos": [0, 1],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[2], "pattern": line_style_list[1],
                 "graph": "g2"},
    "night_mass_cond": {"label": r"$\tau^{\rm cond}_{\rm M,night}$",
                 "label_time": 0.065, "label_pos": [0, 1],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[3], "pattern": line_style_list[3],
                 "graph": "g2"},
    }

time = loadtxt("%s/aging_wc_num_time.txt" % aging_data_dir)

day_num = loadtxt("%s/aging_wc_num_tau_day.txt" % aging_data_dir)
day_num_cond = loadtxt("%s/aging_wc_num_tau_day_cond.txt" % aging_data_dir)
night_num = loadtxt("%s/aging_wc_num_tau_night.txt" % aging_data_dir)
night_num_cond = loadtxt("%s/aging_wc_num_tau_night_cond.txt" % aging_data_dir)

day_mass = loadtxt("%s/aging_wc_mass_tau_day.txt" % aging_data_dir)
day_mass_cond = loadtxt("%s/aging_wc_mass_tau_day_cond.txt" % aging_data_dir)
night_mass = loadtxt("%s/aging_wc_mass_tau_night.txt" % aging_data_dir)
night_mass_cond = loadtxt("%s/aging_wc_mass_tau_night_cond.txt" % aging_data_dir)

day_num = day_num / 3600 # s to hours
day_num_cond = day_num_cond / 3600 # s to hours
night_num = night_num / 3600 # s to hours
night_num_cond = night_num_cond / 3600 # s to hours

day_mass = day_mass / 3600 # s to hours
day_mass_cond = day_mass_cond / 3600 # s to hours
night_mass = night_mass / 3600 # s to hours
night_mass_cond = night_mass_cond / 3600 # s to hours

day_num_smooth = smooth(day_num, window_len = smooth_window_len)
day_num_cond_smooth = smooth(day_num_cond, window_len = smooth_window_len)
night_num_smooth = smooth(night_num, window_len = smooth_window_len)
night_num_cond_smooth = smooth(night_num_cond, window_len = smooth_window_len)

day_mass_smooth = smooth(day_mass, window_len = smooth_window_len)
day_mass_cond_smooth = smooth(day_mass_cond, window_len = smooth_window_len)
night_mass_smooth = smooth(night_mass, window_len = smooth_window_len)
night_mass_cond_smooth = smooth(night_mass_cond, window_len = smooth_window_len)

env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max(time) / 60

for use_color in [True, False]:
    c = canvas.canvas()
    g1 = c.insert(graph.graphxy(
        width = 8.39,
        x = graph.axis.log(min = 0.01,
                              max = 10,
                              title = r"critical supersaturation $S_{\rm c}$ (\%)",
                              painter = grid_painter),
        y = graph.axis.log(min = 1e-2,
                           max = 1e5,
                           title = r"aging timescale $\tau$ (hours)",
                           painter = major_grid_painter)))
    g2 = c.insert(graph.graphxy(
        width = g1.width,
        xpos = g1.xpos + g1.width + grid_h_space,
        x = graph.axis.log(min = 0.01,
                              max = 10,
                              title = r"critical supersaturation $S_{\rm c}$ (\%)",
                              painter = grid_painter),
        y = graph.axis.linkedaxis(g1.axes["y"],
                                  painter = linked_major_grid_painter),
))
#        key = graph.key.key(pos = "br", columns = 2,
#                            keyattrs = [deco.stroked,
#                                        deco.filled([color.rgb.white])])))

    graphs = {"g1": g1, "g2": g2}

    g1.doaxes()
    g2.doaxes()

    x_data = ss_active_axis.edges() * 100

#    for (key, y_data) \
#            in [("day_num", day_num),
#                ("day_num_cond", day_num_cond),
#                ("night_num", night_num),
#                ("night_num_cond", night_num_cond),
#                ("day_mass", day_mass),
#                ("day_mass_cond", day_mass_cond),
#                ("night_mass", night_mass),
#                ("night_mass_cond", night_mass_cond)]:
#        g = graphs[plot_info[key]["graph"]]
#        if use_color:
#            grey_color = color.hsb(plot_info[key]["color"].hsb().color["h"], grey_level, 1)
#            style_attrs = [plot_info[key]["linewidth"],
#                           grey_color]
#        else:
#            grey_color = color.grey(1 - grey_level)
#            style_attrs = [plot_info[key]["linewidth"],
#                           grey_color]
#        plot_data = zip(x_data, y_data)
#        g.plot(
#            graph.data.points(plot_data, x = 1, y = 2),
#            styles = [graph.style.line(lineattrs = style_attrs)])

    for (key, y_data) \
            in [("day_num", day_num),
#                ("day_num_cond", day_num_cond),
                ("night_num", night_num),
                ("night_num_cond", night_num_cond),
                ("day_mass", day_mass),
#                ("day_mass_cond", day_mass_cond),
                ("night_mass", night_mass),
                ("night_mass_cond", night_mass_cond)]:
        if not use_color:
            for supersat in [0.1, 0.3, 0.6, 1.0]:
                i = ss_active_axis.closest_edge(supersat / 100.0)
                timescale = y_data[i]
                print "timescale %s at %f%% supersat = %f h" \
                    % (key, supersat, timescale)
        g = graphs[plot_info[key]["graph"]]
        if use_color:
            style_attrs = [plot_info[key]["linewidth"],
                           plot_info[key]["color"]]
        else:
            style_attrs = [plot_info[key]["linewidth"],
                           plot_info[key]["pattern"]]
        plot_data = zip(x_data, y_data)
        g.plot(
            graph.data.points(plot_data, x = 1, y = 2,
                              title = plot_info[key]["label"]),
            styles = [graph.style.line(lineattrs = style_attrs)])

    for (key, y_data) \
            in [("day_num", day_num),
#                ("day_num_cond", day_num_cond),
                ("night_num", night_num),
                ("night_num_cond", night_num_cond),
                ("day_mass", day_mass),
#                ("day_mass_cond", day_mass_cond),
                ("night_mass", night_mass),
                ("night_mass_cond", night_mass_cond)]:
        g = graphs[plot_info[key]["graph"]]
        plot_data = zip(x_data, y_data)
        label_plot_line_boxed(g, plot_data,
                              plot_info[key]["label_time"],
                              plot_info[key]["label"],
                              plot_info[key]["label_pos"])

    write_text_outside(g1, "day (12:00 LST to 15:00 LST)")
    write_text_outside(g2, "night (18:00 LST to 04:00 LST)")

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    c.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f mm" % unit.tomm(c.bbox().height())
        print "figure width = %.1f mm" % unit.tomm(c.bbox().width())
