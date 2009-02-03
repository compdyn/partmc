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

out_prefix = "figs_aging/aging_aero_time_transfers"

plot_info = {
    "cond_a_f": {"label": r"$\dot{N}^{\rm cond}_{\rm a \to f}$",
                 "label_time": 9.65, "label_pos": [1, 1],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[0], "pattern": line_style_list[1],
                 "graph": "g1"},
    "cond_f_a": {"label": r"$\dot{N}^{\rm cond}_{\rm f \to a}$",
                 "label_time": 10.1, "label_pos": [0, 1],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[2], "pattern": line_style_list[0],
                 "graph": "g1"},
    "coag_loss_a_f": {"label": r"$\dot{N}^{\rm coag}_{\rm a \to f}$",
                 "label_time": 6, "label_pos": [1, 1],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[0], "pattern": line_style_list[1],
                 "graph": "g2"},
    "coag_loss_f_a": {"label": r"$\dot{N}^{\rm coag}_{\rm f \to a}$",
                 "label_time": 13, "label_pos": [0, 1],
                 "linewidth": style.linewidth.Thick,
                 "color": color_list[2], "pattern": line_style_list[0],
                 "graph": "g2"},
    }

time = loadtxt("%s/aging_wc_num_time.txt" % aging_data_dir)
comp_vol = loadtxt("%s/aging_wc_num_comp_vol.txt" % aging_data_dir)

cond_a_f = loadtxt("%s/aging_wc_num_cond_a_f.txt" % aging_data_dir)
cond_f_a = loadtxt("%s/aging_wc_num_cond_f_a.txt" % aging_data_dir)
coag_loss_a_f = loadtxt("%s/aging_wc_num_coag_loss_a_f.txt" % aging_data_dir)
coag_loss_f_a = loadtxt("%s/aging_wc_num_coag_loss_f_a.txt" % aging_data_dir)

cond_a_f_smooth = loadtxt("%s/aging_wc_num_cond_a_f_smooth.txt" % aging_data_dir)
cond_f_a_smooth = loadtxt("%s/aging_wc_num_cond_f_a_smooth.txt" % aging_data_dir)
coag_loss_a_f_smooth = loadtxt("%s/aging_wc_num_coag_loss_a_f_smooth.txt" % aging_data_dir)
coag_loss_f_a_smooth = loadtxt("%s/aging_wc_num_coag_loss_f_a_smooth.txt" % aging_data_dir)

print "level %d = %f%%" % (level_mid, level_mid_value * 100)

cond_a_f_conc = cond_a_f[:,level_mid] / comp_vol[1:] / delta(time) * 1e-6 # m^{-3} s^{-1} to cm^{-3} s^{-1}
cond_f_a_conc = cond_f_a[:,level_mid] / comp_vol[1:] / delta(time) * 1e-6 # m^{-3} s^{-1} to cm^{-3} s^{-1}
coag_loss_a_f_conc = coag_loss_a_f[:,level_mid] / comp_vol[1:] / delta(time) * 1e-6 # m^{-3} s^{-1} to cm^{-3} s^{-1}
coag_loss_f_a_conc = coag_loss_f_a[:,level_mid] / comp_vol[1:] / delta(time) * 1e-6 # m^{-3} s^{-1} to cm^{-3} s^{-1}

cond_a_f_smooth_conc = smooth(cond_a_f_conc, window_len = smooth_window_len)
cond_f_a_smooth_conc = smooth(cond_f_a_conc, window_len = smooth_window_len)
coag_loss_a_f_smooth_conc = smooth(coag_loss_a_f_conc, window_len = smooth_window_len)
coag_loss_f_a_smooth_conc = smooth(coag_loss_f_a_conc, window_len = smooth_window_len)

env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max(time) / 60

for use_color in [True, False]:
    c = canvas.canvas()
    g2 = c.insert(graph.graphxy(
        width = 6.8,
        x = graph.axis.linear(min = 0.,
                              max = max_time_min,
                              parter = graph.axis.parter.linear(tickdists
                                                                = [6 * 60, 3 * 60]),
                              texter = time_of_day(base_time
                                                   = start_time_of_day_min),
                              title = "local standard time (LST) (hours:minutes)",
                              painter = grid_painter),
        y = graph.axis.log(min = 0.002,
                           max = 3,
                           title = r"number conc. rate ($\rm cm^{-3}\,s^{-1}$)",
                           painter = major_grid_painter)))
    g1 = c.insert(graph.graphxy(
        width = 6.8,
        ypos = g2.height + 0.5,
        x = graph.axis.linkedaxis(g2.axes["x"],
                                  painter = graph.axis.painter.linked(gridattrs = [style.linestyle.dotted])),
        y = graph.axis.log(min = 0.002,
                           max = 3,
                           title = r"number conc. rate ($\rm cm^{-3}\,s^{-1}$)",
                           painter = major_grid_painter)))
    graphs = {"g1": g1, "g2": g2}

    for (key, y_data) \
            in [("cond_a_f", cond_a_f_conc),
                ("cond_f_a", cond_f_a_conc),
                ("coag_loss_a_f", coag_loss_a_f_conc),
                ("coag_loss_f_a", coag_loss_f_a_conc)]:
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
        chopped_data = chop_sign_data_iterative(plot_data)
        for signed_data in chopped_data:
            g.plot(
                graph.data.points(signed_data, x = 1, y = 2),
                styles = [graph.style.line(lineattrs = style_attrs)])

    for (key, y_data) \
            in [("cond_a_f", cond_a_f_smooth_conc),
                ("cond_f_a", cond_f_a_smooth_conc),
                ("coag_loss_a_f", coag_loss_a_f_smooth_conc),
                ("coag_loss_f_a", coag_loss_f_a_smooth_conc)]:
        g = graphs[plot_info[key]["graph"]]
        if use_color:
            style_attrs = [plot_info[key]["linewidth"],
                           plot_info[key]["color"]]
        else:
            style_attrs = [plot_info[key]["linewidth"],
                           plot_info[key]["pattern"]]
        plot_data = zip(time[1:] / 60, y_data)
        chopped_data = chop_sign_data_iterative(plot_data)
        for signed_data in chopped_data:
            g.plot(
                graph.data.points(signed_data, x = 1, y = 2),
                styles = [graph.style.line(lineattrs = style_attrs)])

    g1.doaxes()
    g2.doaxes()

    for (key, y_data) \
            in [("cond_a_f", cond_a_f_smooth_conc),
                ("cond_f_a", cond_f_a_smooth_conc),
                ("coag_loss_a_f", coag_loss_a_f_smooth_conc),
                ("coag_loss_f_a", coag_loss_f_a_smooth_conc)]:
        g = graphs[plot_info[key]["graph"]]
        plot_data = zip(time[1:] / 60, y_data)
        label_plot_line_boxed(g, plot_data,
                              plot_info[key]["label_time"] * 60,
                              plot_info[key]["label"],
                              plot_info[key]["label_pos"])

    boxed_text(g1, "condensation", point = [1, 1], anchor_point_rel = [1, 0])
    boxed_text(g2, "coagulation", point = [1, 1], anchor_point_rel = [1, 0])

    write_text_outside(g1, r"critical supersaturation $S = %.1f\%%$" % (level_mid_value * 100))

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    c.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
