#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import time as module_time
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
from fig_helper import *

out_prefix = "movs/aero_2d_all_no_coag_pdfs/aero_2d_all_no_coag"

y_axis_label = r"BC dry mass frac. $w_{{\rm BC},{\rm dry}}\ (\%)$"

def get_plot_data(filename, value_max = None):
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
            for j in range(len(show_particles)):
                if particles.id[i] == show_particles[j]["id"]:
                    show_coords[j] = [diameter[i], comp_frac[i]]

    value = num_den_array / num_den_array.sum() \
            / x_axis.grid_size(0) / (y_axis.grid_size(0) / 100.0)
    if value_max == None:
        value_max = value.max()
    if value_max > 0.0:
        value = value / value_max
    value = value.clip(0.0, 1.0)

    rects = pmc_histogram_2d_multi([value],
                                    x_axis, y_axis)
    return (rects, show_coords, env_state)

time_filename_list = get_time_filename_list(netcdf_dir_nc, netcdf_pattern_nc)
count = 0
start_time = module_time.clock()
for [time, filename, output_key] in time_filename_list:
    g = make_1x1_graph_grid(y_axis_label)
    (rects, show_coords, env_state) = get_plot_data(filename, max_val)
    g.plot(graph.data.points(rects,
                             xmin = 1, xmax = 2, ymin = 3, ymax = 4,
                             color = 5),
           styles = [hsb_rect(rainbow_palette)])

    write_time(g, env_state, with_hours = False)

    g.dolayout()
    for axisname in ["x", "y"]:
        for t in g.axes[axisname].data.ticks:
            if t.ticklevel is not None:
                g.stroke(g.axes[axisname].positioner.vgridpath(t.temp_v),
                         [style.linestyle.dotted])
    g.dodata()
    g.doaxes()

    add_canvas_color_bar(g,
                         min = 0.0,
                         max = max_val,
                         xpos = g.xpos + g.width + grid_h_space,
                         ybottom = g.ypos,
                         ytop = g.ypos + g.height,
                         title = r"normalized number conc. $\hat{n}_{\rm BC,dry}(D,w)$",
                         palette = rainbow_palette)

    out_filename = "%s_%s.pdf" % (out_prefix, output_key)
    g.writePDFfile(out_filename)

    count += 1
    total_count = len(time_filename_list)
    elapsed_time = module_time.clock() - start_time
    total_time = float(total_count) / float(count) * elapsed_time
    remaining_time = total_time - elapsed_time
    print "all, no coag: %d of %d - %.1f s elapsed, %.1f s remaining" \
          % (count, total_count, elapsed_time, remaining_time)
