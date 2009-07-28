#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import random as py_random
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append(".")
from fig_helper import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *

out_prefix = "figs_aging/aging_aero_2d_age_disc"

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

def percentile(data, p):
    # p in [0,1]
    # data must be sorted
    i = int(floor((len(data) - 1) * p + 0.5))
    return data[i]

def grid_plot_data(x, y, z, x_axis, y_axis, hi_x_axis, hi_y_axis):
    x_bin = x_axis.find(x)
    y_bin = y_axis.find(y)
    values = [[[] for j in range(y_axis.n_bin)]
              for i in range(x_axis.n_bin)]
    for i in range(len(x)):
        if x_axis.valid_bin(x_bin[i]) and y_axis.valid_bin(y_bin[i]):
            values[x_bin[i]][y_bin[i]].append(z[i])
    for x_bin in range(x_axis.n_bin):
        for y_bin in range(y_axis.n_bin):
            values[x_bin][y_bin].sort()
    grid = numpy.zeros([hi_x_axis.n_bin, hi_y_axis.n_bin])
    mask = numpy.zeros([hi_x_axis.n_bin, hi_y_axis.n_bin], int)
    for x_bin in range(x_axis.n_bin):
        for y_bin in range(y_axis.n_bin):
            if len(values[x_bin][y_bin]) > 0:
                subs = [(0,0),(0,1),(1,0),(1,1)]
                py_random.shuffle(subs)
                (sub_min, sub_max, sub_low, sub_high) = subs
                val_min = min(values[x_bin][y_bin])
                val_max = max(values[x_bin][y_bin])
                val_low = percentile(values[x_bin][y_bin], 0.3333)
                val_high = percentile(values[x_bin][y_bin], 0.6666)
                for (sub, val) in [(sub_min, val_min),
                                   (sub_max, val_max),
                                   (sub_low, val_low),
                                   (sub_high, val_high)]:
                    sub_i, sub_j = sub
                    i = x_bin * 2 + sub_i
                    j = y_bin * 2 + sub_j
                    grid[i,j] = val
                    mask[i,j] = 1
    plot_data = pmc_histogram_2d(grid, hi_x_axis, hi_y_axis, mask = mask)
    return plot_data

time_filename_list = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)

time = 24 * 3600.0
filename = file_filename_at_time(time_filename_list, time)

ncf = NetCDFFile(filename)
particles = aero_particle_array_t(ncf)
env_state = env_state_t(ncf)
ncf.close()

diameter = particles.dry_diameter() * 1e6
critical_ss = (particles.kappa_rh(env_state, const) - 1.0) * 100.0
comp_frac = particles.mass(include = ["BC"]) \
            / particles.mass(exclude = ["H2O"]) * 100
# hack to avoid landing just around the integer boundaries
comp_frac *= (1.0 + 1e-12)
age = (env_state.elapsed_time - particles.least_create_time) / 3600.0
max_age = 24.0
scaled_age = age / max_age

diameter_axis = pmc_log_axis(min = diameter_axis_min,
                             max = diameter_axis_max,
                             n_bin = num_diameter_bins)
hi_diameter_axis = pmc_log_axis(min = diameter_axis_min,
                             max = diameter_axis_max,
                             n_bin = num_diameter_bins * 2)
bc_axis = pmc_linear_axis(min = bc_axis_min, max = bc_axis_max,
                          n_bin = num_bc_bins)
hi_bc_axis = pmc_linear_axis(min = bc_axis_min, max = bc_axis_max,
                          n_bin = num_bc_bins * 2)
ss_axis = pmc_log_axis(min = ss_axis_min, max = ss_axis_max,
                       n_bin = num_ss_bins)
hi_ss_axis = pmc_log_axis(min = ss_axis_min, max = ss_axis_max,
                       n_bin = num_ss_bins * 2)

ss_plot_data = grid_plot_data(diameter, critical_ss, scaled_age,
                              diameter_axis, ss_axis,
                              hi_diameter_axis, hi_ss_axis)
bc_plot_data = grid_plot_data(diameter, comp_frac, scaled_age,
                              diameter_axis, bc_axis,
                              hi_diameter_axis, hi_bc_axis)

for color in [True, False]:
    c = canvas.canvas()
    g1 = c.insert(graph.graphxy(
        width = grid_graph_width,
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
                           title = diameter_axis_label),
        y = graph.axis.log(min = ss_axis_min,
                              max = ss_axis_max,
                              title = r"critical supersaturation $S_{\rm c}$ (\%)")))
    g2 = c.insert(graph.graphxy(
            width = grid_graph_width,
            ypos = g1.height + grid_v_space,
            x = graph.axis.linkedaxis(g1.axes["x"]),
            y = graph.axis.linear(min = bc_axis_min,
                                  max = bc_axis_max,
                                  title = r"BC dry mass frac. $w_{{\rm BC},{\rm dry}}\ (\%)$",
                                  density = 1.2)))

    if color:
        palette = rainbow_palette
    else:
        palette = nonlinear_gray_palette

    for g in [g1, g2]:
        g.dolayout()

        for axisname in ["x", "y"]:
            for t in g.axes[axisname].data.ticks:
                if t.ticklevel is not None:
                    g.stroke(g.axes[axisname].positioner.vgridpath(t.temp_v),
                             [style.linestyle.dotted])
                    
    for (plot_data, g) in [(ss_plot_data, g1),
                           (bc_plot_data, g2)]:
        g.plot(graph.data.points(plot_data,
                                 xmin = 1, xmax = 2, ymin = 3, ymax = 4,
                                 color = 5),
               styles = [hsb_rect(palette)])

    for g in [g1, g2]:
        write_time(g, env_state)

        g.dodata()
        g.doaxes()

        add_canvas_color_bar(
            g,
            min = 0.0,
            max = max_age,
            xpos = g.xpos + g.width + grid_h_space,
            ybottom = g.ypos,
            ytop = g.ypos + g.height,
            title = r"age (hours)",
            palette = palette,
            parter = graph.axis.parter.linear(tickdists = [6, 3]))


    if color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    c.writePDFfile(out_filename)
    if not color:
        print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
