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

out_prefix = "figs_aging/aging_aero_time_totals"

plot_info = {
    "total": {"label": r"$N_{\rm tot}$", "label_time": 8, "label_pos": [1, 1],
              "linewidth": style.linewidth.Thick,
              "color": color_list[0], "pattern": line_style_list[0]},
    "aged":  {"label": r"$N_{\rm a}$", "label_time": 15, "label_pos": [1, 0],
              "linewidth": style.linewidth.Thick,
              "color": color_list[1], "pattern": line_style_list[1]},
    "fresh": {"label": r"$N_{\rm f}$", "label_time": 15, "label_pos": [0, 1],
              "linewidth": style.linewidth.Thick,
              "color": color_list[2], "pattern": line_style_list[2]},
    }

time = loadtxt("%s/aging_wc_num_time.txt" % aging_data_dir)
comp_vol = loadtxt("%s/aging_wc_num_comp_vol.txt" % aging_data_dir)
num_a = loadtxt("%s/aging_wc_num_a.txt" % aging_data_dir)
num_f = loadtxt("%s/aging_wc_num_f.txt" % aging_data_dir)
num_t = num_a + num_f

print "level %d = %f%%" % (level_mid, level_mid_value * 100)
num_a_conc = num_a[:,level_mid] / comp_vol * 1e-6 # m^{-3} to cm^{-3}
num_f_conc = num_f[:,level_mid] / comp_vol * 1e-6 # m^{-3} to cm^{-3}
num_t_conc = num_t[:,level_mid] / comp_vol * 1e-6 # m^{-3} to cm^{-3}

env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max(time) / 60

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
                              max = 1.5e4,
                              #density = 1.5, # hack to increase number of ticks
                              title = r"number conc. ($\rm cm^{-3}$)",
                              painter = grid_painter))

    g.doaxes()

    for (key, y_data) in [("aged", num_a_conc),
                          ("fresh", num_f_conc),
                          ("total", num_t_conc)]:
        if use_color:
            style_attrs = [plot_info[key]["linewidth"],
                           plot_info[key]["color"]]
        else:
            style_attrs = [plot_info[key]["linewidth"],
                           plot_info[key]["pattern"]]
        plot_data = zip(time / 60, y_data)
        g.plot(
            graph.data.points(plot_data, x = 1, y = 2),
            styles = [graph.style.line(lineattrs = style_attrs)])
        label_plot_line_boxed(g, plot_data,
                              plot_info[key]["label_time"] * 60,
                              plot_info[key]["label"],
                              plot_info[key]["label_pos"])

    write_text_outside(g, r"critical supersaturation $S_{\rm c} = %.1f\%%$" % (level_mid_value * 100))

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    g.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
