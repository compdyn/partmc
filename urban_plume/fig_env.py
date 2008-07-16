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

out_filename = "figs/env.pdf"

netcdf_dir = "out"
netcdf_pattern = r"urban_plume_state_0001_([0-9]{8})\.nc"

env_state_history = read_history(env_state_t, netcdf_dir, netcdf_pattern)
start_time_of_day_min = env_state_history[0][1].start_time_of_day / 60
max_time_min = max([time for [time, env_state] in env_state_history]) / 60

g = graph.graphxy(
    width = 10,
    height = 4,
    x = graph.axis.linear(min = 0.,
                          max = max_time_min,
			  title = "local standard time (hours:minutes)",
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time
                                               = start_time_of_day_min),
			  painter = grid_painter),
    y = graph.axis.linear(min = 285,
                          max = 300,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [3, 1.5]),
                          title = "temperature (K)",
                          painter = grid_painter),
    y2 = graph.axis.linear(min = 50,
                           max = 100,
                           parter = graph.axis.parter.linear(tickdists
                                                             = [10, 5]),
                           title = "relative humidity (1)",
                           texter = graph.axis.texter.decimal(suffix = r"\%")),
    y4 = graph.axis.linear(min = 0,
                           max = 500,
                           parter = graph.axis.parter.linear(tickdists
                                                             = [100, 50]),
                           title = "mixing height (m)"))

temp_plot_data = []
rh_plot_data = []
height_plot_data = []
for [time, env_state] in env_state_history:
    temp_plot_data.append([time / 60, env_state.temperature])
    rh_plot_data.append([time / 60, env_state.relative_humidity * 100])
    height_plot_data.append([time / 60, env_state.height])

g.plot(graph.data.points(temp_plot_data, x = 1, y = 2),
       styles = [graph.style.line(lineattrs = [line_style_list[0],
                                               style.linewidth.THick])])
g.plot(graph.data.points(rh_plot_data, x = 1, y2 = 2),
       styles = [graph.style.line(lineattrs = [line_style_list[1],
                                               style.linewidth.THick])])
g.plot(graph.data.points(height_plot_data, x = 1, y4 = 2),
       styles = [graph.style.line(lineattrs = [line_style_list[2],
                                               style.linewidth.THick])])

label_plot_line(g, temp_plot_data, 8 * 60.0, "temperature", [0, 1],
                1 * unit.v_mm)
label_plot_line(g, rh_plot_data, 8 * 60.0, "relative humidity", [0, 0],
                1 * unit.v_mm, yaxis = g.axes["y2"])
label_plot_line(g, height_plot_data, 16 * 60.0, "mixing height", [0, 1],
                1 * unit.v_mm, yaxis = g.axes["y4"])

g.writePDFfile(out_filename)
