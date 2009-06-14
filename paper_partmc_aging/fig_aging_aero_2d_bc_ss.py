#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append(".")
from fig_helper import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *

out_prefix = "figs_aging/aging_aero_2d_bc_ss"

y_axis_label = r"BC dry mass frac. $w_{{\rm BC},{\rm dry}}\ (\%)$"

ss_axis_min = 1e-2
ss_axis_max = 1e2
num_ss_bins = 40

bc_max_val = 4.0
ss_max_val = 100.0

const = load_constants("../src/constants.f90")

def get_plot_data_bc(filename, value_max = None):
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()

    diameter = particles.dry_diameter() * 1e6
    comp_frac = particles.mass(include = ["BC"]) \
                / particles.mass(exclude = ["H2O"]) * 100

    x_axis = pmc_log_axis(min = diameter_axis_min, max = diameter_axis_max,
                          n_bin = num_diameter_bins)
    y_axis = pmc_linear_axis(min = bc_axis_min, max = bc_axis_max,
                             n_bin = num_bc_bins)
    x_bin = x_axis.find(diameter)
    # hack to avoid landing just around the integer boundaries
    comp_frac *= (1.0 + 1e-12)
    y_bin = y_axis.find(comp_frac)

    num_den_array = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
    show_coords = [[] for p in show_particles]
    for i in range(particles.n_particles):
        if x_axis.valid_bin(x_bin[i]) and y_axis.valid_bin(y_bin[i]):
            scale = particles.comp_vol[i] * x_axis.grid_size(x_bin[i]) \
                    * (y_axis.grid_size(y_bin[i]) / 100)
            num_den_array[x_bin[i], y_bin[i]] += 1.0 / scale

    value = num_den_array / num_den_array.sum() \
            / x_axis.grid_size(0) / (y_axis.grid_size(0) / 100.0)
    if value_max == None:
        value_max = value.max()
    if value_max > 0.0:
        value = value / value_max
    value = value.clip(0.0, 1.0)

    rects = pmc_histogram_2d_multi([value],
                                    x_axis, y_axis)
    return (rects, env_state)

def get_plot_data_ss(filename, value_max = None):
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()

    diameter = particles.dry_diameter() * 1e6
    critical_ss = (particles.kappa_rh(env_state, const) - 1.0) * 100.0

    x_axis = pmc_log_axis(min = diameter_axis_min, max = diameter_axis_max,
                          n_bin = num_diameter_bins)
    y_axis = pmc_log_axis(min = ss_axis_min, max = ss_axis_max,
                          n_bin = num_ss_bins)
    x_bin = x_axis.find(diameter)
    # hack to avoid landing just around the integer boundaries
    #comp_frac *= (1.0 + 1e-12)
    y_bin = y_axis.find(critical_ss)

    num_den_array = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
    show_coords = [[] for p in show_particles]
    for i in range(particles.n_particles):
        if x_axis.valid_bin(x_bin[i]) and y_axis.valid_bin(y_bin[i]):
            scale = particles.comp_vol[i] * x_axis.grid_size(x_bin[i]) \
                    * (y_axis.grid_size(y_bin[i]) / 100)
            num_den_array[x_bin[i], y_bin[i]] += 1.0 / scale

    value = num_den_array / num_den_array.sum() \
            / x_axis.grid_size(0) / (y_axis.grid_size(0) / 100.0)
    if value_max == None:
        value_max = value.max()
    if value_max > 0.0:
        value = value / value_max
    value = value.clip(0.0, 1.0)

    rects = pmc_histogram_2d_multi([value],
                                    x_axis, y_axis)
    return (rects, env_state)

time_filename_list = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
for color in [True, False]:
    c = canvas.canvas()
    g21 = c.insert(graph.graphxy(
        width = grid_graph_width,
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
                           title = diameter_axis_label),
        y = graph.axis.log(min = ss_axis_min,
                              max = ss_axis_max,
                              title = r"critical supersaturation $S_{\rm c}$ (\%)")))
    g11 = c.insert(graph.graphxy(
            width = grid_graph_width,
            ypos = g21.height + grid_v_space,
            x = graph.axis.linkedaxis(g21.axes["x"]),
            y = graph.axis.linear(min = bc_axis_min,
                                  max = bc_axis_max,
                                  title = r"BC dry mass frac. $w_{{\rm BC},{\rm dry}}\ (\%)$",
                                  density = 1.2)))
    g22 = c.insert(graph.graphxy(
        width = grid_graph_width,
        xpos = g21.width + grid_h_space,
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
                           title = diameter_axis_label),
        y = graph.axis.linkedaxis(g21.axes["y"])))
    g12 = c.insert(graph.graphxy(
        width = grid_graph_width,
        xpos = g11.width + grid_h_space,
        ypos = g22.height + grid_v_space,
        x = graph.axis.linkedaxis(g22.axes["x"]),
        y = graph.axis.linkedaxis(g11.axes["y"])))

    graphs = {"c": c,
              "g11" : g11,
              "g12" : g12,
              "g21" : g21,
              "g22" : g22}

    plot_info_list = [
        {"graph": "g11", "y_axis": "bc", "time_hour": 1},
        {"graph": "g12", "y_axis": "bc", "time_hour": 24},
        {"graph": "g21", "y_axis": "ss", "time_hour": 1},
        {"graph": "g22", "y_axis": "ss", "time_hour": 24}]

    for plot_info in plot_info_list:
        time = plot_info["time_hour"] * 3600.0
        filename = file_filename_at_time(time_filename_list, time)
        if plot_info["y_axis"] == "bc":
            (rects, env_state) = get_plot_data_bc(filename, bc_max_val)
        else:
            (rects, env_state) = get_plot_data_ss(filename, ss_max_val)
        g = graphs[plot_info["graph"]]
        if color:
            palette = rainbow_palette
        else:
            palette = gray_palette
        g.plot(graph.data.points(rects,
                                 xmin = 1, xmax = 2, ymin = 3, ymax = 4,
                                 color = 5),
               styles = [hsb_rect(palette)])

        write_time(g, env_state)

        g.dolayout()
        for axisname in ["x", "y"]:
            for t in g.axes[axisname].data.ticks:
                if t.ticklevel is not None:
                    g.stroke(g.axes[axisname].positioner.vgridpath(t.temp_v),
                             [style.linestyle.dotted])
        g.dodata()
        g.doaxes()

    c = graphs["c"]
    add_canvas_color_bar(c,
                         min = 0.0,
                         max = bc_max_val,
                         xpos = graphs["g12"].xpos + graphs["g12"].width + grid_h_space,
                         ybottom = graphs["g12"].ypos,
                         ytop = graphs["g12"].ypos + graphs["g12"].height,
                         title = r"normalized number conc.",
                         palette = palette)
    add_canvas_color_bar(c,
                         min = 0.0,
                         max = ss_max_val,
                         xpos = graphs["g22"].xpos + graphs["g22"].width + grid_h_space,
                         ybottom = graphs["g22"].ypos,
                         ytop = graphs["g22"].ypos + graphs["g22"].height,
                         title = r"normalized number conc.",
                         palette = palette)

    if color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    c.writePDFfile(out_filename)
    if not color:
        print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
