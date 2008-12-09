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
sys.path.append(".")
from fig_helper import *

out_prefix = "figs/env"

env_state_history = read_history(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
start_time_of_day_min = env_state_history[0][1].start_time_of_day / 60
max_time_min = max([time for [time, env_state] in env_state_history]) / 60

for use_color in [True, False]:
    g = graph.graphxy(
        width = 10,
        height = 4,
        x = graph.axis.linear(min = 0.,
                              max = max_time_min,
                              title = "local standard time (LST) (hours:minutes)",
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
                               title = "relative humidity ($\%$)"),
        y4 = graph.axis.linear(min = 0,
                               max = 500,
                               parter = graph.axis.parter.linear(tickdists
                                                                 = [100, 50]),
                               title = "mixing height (m)"))

    g.doaxes()

    temp_plot_data = []
    rh_plot_data = []
    height_plot_data = []
    for [time, env_state] in env_state_history:
        temp_plot_data.append([time / 60, env_state.temperature])
        rh_plot_data.append([time / 60, env_state.relative_humidity * 100])
        height_plot_data.append([time / 60, env_state.height])

    if use_color:
        use_line_style_list = color_list
    else:
        use_line_style_list = line_style_list
    g.plot(graph.data.points(temp_plot_data, x = 1, y = 2),
           styles = [graph.style.line(lineattrs = [use_line_style_list[0],
                                                   style.linewidth.THick])])
    g.plot(graph.data.points(rh_plot_data, x = 1, y2 = 2),
           styles = [graph.style.line(lineattrs = [use_line_style_list[1],
                                                   style.linewidth.THick])])
    g.plot(graph.data.points(height_plot_data, x = 1, y4 = 2),
           styles = [graph.style.line(lineattrs = [use_line_style_list[2],
                                                   style.linewidth.THick])])

    label_plot_line_boxed(g, temp_plot_data, 9.7 * 60.0, "temperature", [0, 1])
    label_plot_line_boxed(g, rh_plot_data, 4 * 60.0, "relative humidity", [0, 1],
                    yaxis = g.axes["y2"])
    label_plot_line_boxed(g, height_plot_data, 22 * 60.0, "mixing height", [1, 0],
                    yaxis = g.axes["y4"])

    if not use_color:
        print "init RH = %g%%" % rh_plot_data[0][1]
        print "min RH = %g%%" % min([v for [t,v] in rh_plot_data])
        print "max RH in night = %g%%" % max([v for [t,v] in rh_plot_data
                                              if t > 12 * 60])
        t_85 = min([t for [t,v] in rh_plot_data
                    if v < 85])
        i = find_nearest_time(rh_plot_data, t_85)
        print "RH at %s = %g%%" \
              % (time_of_day_string(rh_plot_data[i][0] * 60
                                    + env_state_history[0][1].start_time_of_day),
                 rh_plot_data[i][1])

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    g.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
