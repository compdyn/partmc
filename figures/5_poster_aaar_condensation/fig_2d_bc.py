#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../../tool")
from pmc_data_nc import *
from pmc_pyx import *
from config import *

ss_axis_min = 1e-2
ss_axis_max = 1e2
num_ss_bins = 40

const = load_constants("../../src/constants.f90")

out_prefix = "figs/2d_bc"

def get_plot_data_bc(filename, value_min = None, value_max = None):
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
    for i in range(particles.n_particles):
        if x_axis.valid_bin(x_bin[i]) and y_axis.valid_bin(y_bin[i]):
            scale = particles.comp_vol[i] * x_axis.grid_size(x_bin[i]) \
                    * (y_axis.grid_size(y_bin[i]) / 100)
            num_den_array[x_bin[i], y_bin[i]] += 1.0 / scale

    value = num_den_array / 1e6
    #value = num_den_array / num_den_array.sum() \
    #        / x_axis.grid_size(0) / (y_axis.grid_size(0) / 100.0)
    if value_max == None:
        value_max = value.max()
    if value_min == None:
        maxed_value = where(value > 0.0, value, value_max)
        value_min = maxed_value.min()
    if value_max > 0.0:
        value = (log(value) - log(value_min)) \
                / (log(value_max) - log(value_min))
    value = value.clip(0.0, 1.0)

    rects = pmc_histogram_2d_multi([value],
                                    x_axis, y_axis)
    return (rects, env_state, value_min, value_max)

for [i_run, netcdf_pattern] in netcdf_indexed_patterns:
    out_filename = "%s_%d.pdf" % (out_prefix, i_run)
    print out_filename

    g = graph.graphxy(
        width = grid_graph_width,
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
                           title = diameter_axis_label),
        y = graph.axis.linear(min = bc_axis_min,
                              max = bc_axis_max,
                              title = r"BC dry mass frac. $w_{{\rm BC},{\rm dry}}\ (\%)$",
                              density = 1.2))

    filename_list = get_filename_list(netcdf_dir, netcdf_pattern)
    filename = filename_list[0]
    (rects, env_state, bc_min_val, bc_max_val) \
        = get_plot_data_bc(filename)

    palette = rainbow_palette

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

    add_canvas_color_bar(g,
                         min = bc_min_val,
                         max = bc_max_val,
                         log_scale = True,
                         xpos = g.xpos + g.width + grid_h_space,
                         ybottom = g.ypos,
                         ytop = g.ypos + g.height,
                         title = r"number conc. $(\rm cm^{-3})$",
                         palette = palette)

    g.writePDFfile(out_filename)
    print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
    print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
