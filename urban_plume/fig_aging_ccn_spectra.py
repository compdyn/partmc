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

plot_info_list = [
    {"time_hour": 1, "label_x": 1.2, "label_pos": [0, 0],
     "linewidth": style.linewidth.Thick,
     "color": color_list[0], "pattern": line_style_list[0]},
    {"time_hour": 24, "label_x": 0.18, "label_pos": [1, 1],
     "linewidth": style.linewidth.Thick,
     "color": color_list[2], "pattern": line_style_list[1]},
    ]

out_prefix = "figs_aging/aging_ccn_spectra"

const = load_constants("../src/constants.f90")

time_filename_list = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max([time for [time, filename, key] in time_filename_list]) / 60

def get_plot_data(filename):
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()

    critical_ss = (particles.kappa_rh(env_state, const) - 1.0) * 100.0
    critical_ss.sort()
    plot_data = cumulative_hi_res(
        critical_ss,
        ones_like(critical_ss) / critical_ss.size * 100.0,
        final = 100.0,
        min_x_factor = 1.1, min_y_step = 0.3)

    return (plot_data, env_state)

for use_color in [True, False]:
    g = graph.graphxy(
        width = 6.8,
        x = graph.axis.log(min = 0.01,
                           max = 10,
                           title = r"critical supersaturation $S$ (\%)",
                           painter = grid_painter),
        y = graph.axis.linear(min = 0.0,
                              max = 100.0,
                              title = r"cumulative number fraction (\%)",
                              painter = grid_painter))

    plot_data_list = []
    for plot_info in plot_info_list:
        time = plot_info["time_hour"] * 3600.0
        filename = file_filename_at_time(time_filename_list, time)
        (plot_data, env_state) = get_plot_data(filename)
        plot_data_list.append(plot_data)
        if use_color:
            style_attrs = [plot_info["linewidth"],
                           plot_info["color"]]
        else:
            style_attrs = [plot_info["linewidth"],
                           plot_info["pattern"]]
        g.plot(
            graph.data.points(plot_data, x = 1, y = 2),
            styles = [graph.style.line(lineattrs = style_attrs)])
    
    g.dodata()
    g.doaxes()

    for [plot_info, plot_data] in zip(plot_info_list, plot_data_list):
        if plot_info["time_hour"] == 1:
            label = "%d hour" % plot_info["time_hour"]
        else:
            label = "%d hours" % plot_info["time_hour"]
        label_plot_line_boxed(g, plot_data, plot_info["label_x"],
                              label, plot_info["label_pos"])

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    g.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
