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

out_prefix = "figs/aero_water_dist"

time_filename_list = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
env_state_history = read_history(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time = max([time for [time, filename, key] in time_filename_list])
max_time_min = max_time / 60

for use_color in [True, False]:
    g = graph.graphxy(
        width = 6.7,
        x = graph.axis.linear(min = 0.0,
                              max = max_time_min,
                              parter = graph.axis.parter.linear(tickdists
                                                                = [6 * 60, 3 * 60]),
                              texter = time_of_day(base_time
                                                   = start_time_of_day_min),
                              title = "local standard time (LST) (hours:minutes)",
                              painter = grid_painter),
        y = graph.axis.linear(min = 0.,
                              max = 100,
                              parter = graph.axis.parter.linear(tickdists
                                                                = [20, 10]),
                              title = r"water mass frac. $w_{\rm H_2O,all}$ ($\%$)",
                              painter = grid_painter),
        y2 = graph.axis.linear(min = 0,
                               max = 100,
                               parter = graph.axis.parter.linear(tickdists
                                                                 = [20, 10]),
                               title = "relative humidity ($\%$)"))

    x_axis = pmc_linear_axis(min = 0, max = max_time_min, n_bin = 500)
    y_axis = pmc_linear_axis(min = 0, max = 100, n_bin = 100)
    num_den_array = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
    num_times_array = numpy.zeros([x_axis.n_bin], dtype=int)

    for [time, filename, key] in time_filename_list:
        ncf = NetCDFFile(filename)
        particles = aero_particle_array_t(ncf)
        ncf.close()

        water_frac = particles.mass(include = ["H2O"]) \
                     / particles.mass() * 100
        x_bin = x_axis.find(array([time / 60.0]))[0]
        y_bin = y_axis.find(water_frac)

        num_times_array[x_bin] += 1
        for i in range(particles.n_particles):
            scale = particles.comp_vol[i] * (y_axis.grid_size(y_bin[i]) / 100)
            num_den_array[x_bin, y_bin[i]] += 1.0 / scale * 1e-6 # m^{-3} to cm^{-3}

    for i in range(x_axis.n_bin):
        num_den_array[i,:] /= num_times_array[i]

    value = num_den_array
    value_max = value.max() / 10.0
    if value_max != None:
        value = value / value_max
    value = value.clip(0.0, 1.0)

    rects = pmc_histogram_2d_multi([value], x_axis, y_axis)

    if use_color:
        palette = rainbow_palette
    else:
        palette = gray_palette
    g.plot(graph.data.points(rects,
                             xmin = 1, xmax = 2, ymin = 3, ymax = 4,
                             color = 5),
           styles = [hsb_rect(palette)])

    rh_plot_data = []
    for [time, env_state] in env_state_history:
        rh_plot_data.append([time / 60, env_state.relative_humidity * 100])
    if use_color:
        rh_attr = style.linestyle.solid
    else:
        rh_attr = style.linestyle.dashed
    g.plot(graph.data.points(rh_plot_data, x = 1, y2 = 2),
           styles = [graph.style.line(lineattrs = [color.gray.white,
                                                   style.linewidth.THICk]),
                     graph.style.line(lineattrs = [color.gray.black,
                                                   rh_attr,
                                                   style.linewidth.THick])])

    g.dolayout()
    for axisname in ["x", "y"]:
        for t in g.axes[axisname].data.ticks:
            if t.ticklevel is not None:
                g.stroke(g.axes[axisname].positioner.vgridpath(t.temp_v),
                         [style.linestyle.dotted])
    g.dodata()
    g.doaxes()

    label_plot_line_boxed(g, rh_plot_data, 11.5 * 60.0, "relative humidity", [1, 1],
                          yaxis = g.axes["y2"])

    add_color_bar(g,
                  min = 0.0,
                  max = value_max,
                  title = r"number conc. ($\rm cm^{-3}$)",
                  palette = palette,
                  bar_x_offset = 2.5)

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    g.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
