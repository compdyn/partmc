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

netcdf_dir = "out"
netcdf_pattern = r"urban_plume_state_0001_([0-9]{8})\.nc"

env_state_history = read_history(env_state_t, netcdf_dir, netcdf_pattern)
start_time_of_day_min = env_state_history[0][1].start_time_of_day / 60
max_time_min = max([time for [time, env_state] in env_state_history]) / 60

g = graph.graphxy(
    width = 10,
    x = graph.axis.linear(min = 0.0,
                          max = max_time_min,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time
                                               = start_time_of_day_min),
                          title = "time (LST)",
			  painter = grid_painter),
    y = graph.axis.linear(title = "temperature (K)"),
    y2 = graph.axis.linear(title = r"relative humidity (1)",
                           texter = graph.axis.texter.decimal(suffix = r"\%")),
    y4 = graph.axis.linear(title = "mixing height (m)"),
    key = graph.key.key(pos = "tr"))

temp_plot_data = []
rh_plot_data = []
height_plot_data = []
for [time, env_state] in env_state_history:
    temp_plot_data.append([time / 60, env_state.temperature])
    rh_plot_data.append([time / 60, env_state.relative_humidity * 100])
    height_plot_data.append([time / 60, env_state.height])

g.plot(graph.data.points(temp_plot_data, x = 1, y = 2,
                         title = "temperature"),
       styles = [graph.style.line(lineattrs = [color_list[0],
                                               style.linewidth.THick])])
g.plot(graph.data.points(rh_plot_data, x = 1, y2 = 2,
                         title = "relative humidity"),
       styles = [graph.style.line(lineattrs = [color_list[1],
                                               style.linewidth.THick])])
g.plot(graph.data.points(height_plot_data, x = 1, y4 = 2,
                         title = "mixing height"),
       styles = [graph.style.line(lineattrs = [color_list[2],
                                               style.linewidth.THick])])

g.writePDFfile("out/env.pdf")
