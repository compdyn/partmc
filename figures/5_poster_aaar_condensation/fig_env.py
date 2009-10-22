#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../../tool")
from pmc_data_nc import *
from pmc_pyx import *
from config import *

out_prefix = "figs/env"


for [i_run, netcdf_pattern] in netcdf_indexed_patterns:
    out_filename = "%s_%d.pdf" % (out_prefix, i_run)
    print out_filename

    env_state_history = read_history(env_state_t, netcdf_dir, netcdf_pattern)
    start_time_of_day_min = env_state_history[0][1].start_time_of_day / 60
    min_time_min = min([time for [time, env_state] in env_state_history]) / 60
    max_time_min = max([time for [time, env_state] in env_state_history]) / 60

    g = graph.graphxy(
        width = graph_width,
        x = graph.axis.linear(min = 0,
                              max = max_time_min - min_time_min,
                              title = "time (min)",
                              painter = grid_painter),
        y = graph.axis.linear(min = 285,
                              max = 290,
                              title = "temperature (K)",
                              painter = grid_painter),
        y2 = graph.axis.linear(min = 0,
                               max = 0.2,
                               title = "supersaturation ($\%$)",
                               parter = graph.axis.parter.linear(tickdists
                                                                = [0.04, 0.02])))

    g.doaxes()

    temp_plot_data = []
    rh_plot_data = []
    for [time, env_state] in env_state_history:
        temp_plot_data.append([time / 60 - min_time_min, env_state.temperature])
        rh_plot_data.append([time / 60 - min_time_min, (env_state.relative_humidity - 1) * 100])

    use_line_style_list = color_list
    g.plot(graph.data.points(temp_plot_data, x = 1, y = 2),
           styles = [graph.style.line(lineattrs = [use_line_style_list[0],
                                                   style.linewidth.Thick])])
    g.plot(graph.data.points(rh_plot_data, x = 1, y2 = 2),
           styles = [graph.style.line(lineattrs = [use_line_style_list[1],
                                                   style.linewidth.Thick])])

    label_plot_line_boxed(g, temp_plot_data, 4.5, "temperature", [0, 1])
    label_plot_line_boxed(g, rh_plot_data, 5.5, "supersaturation", [1, 0],
                          yaxis = g.axes["y2"])

    print "init supersat = %g%%" % rh_plot_data[0][1]
    print "max supersat = %.10g%%" % max([v for [t,v] in rh_plot_data])

    #write_time(g, env_state_history[0][1])

    g.writePDFfile(out_filename)
    print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
    print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
