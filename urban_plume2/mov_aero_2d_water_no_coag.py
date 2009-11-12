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

out_prefix = "movs/aero_2d_water_no_coag_pdfs/aero_2d_water_no_coag"

y_axis_label = r"BC dry mass frac. $w_{{\rm BC},{\rm dry}}\ (\%)$"

max_val = 50
min_palette_index = 0

def get_plot_data(filename, value_max = None, print_info = True):
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()

    diameter = particles.dry_diameter() * 1e6
    comp_frac = particles.mass(include = ["BC"]) \
                 / particles.mass(exclude = ["H2O"]) * 100
    water_frac = particles.mass(include = ["H2O"]) \
                 / particles.mass() * 100

    x_axis = pmc_log_axis(min = diameter_axis_min, max = diameter_axis_max,
                          n_bin = num_diameter_bins)
    y_axis = pmc_linear_axis(min = bc_axis_min, max = bc_axis_max,
                             n_bin = num_bc_bins)
    x_bin = x_axis.find(diameter)
    # hack to avoid landing just around the integer boundaries
    comp_frac *= (1.0 + 1e-12)
    y_bin = y_axis.find(comp_frac)

    water_frac_array = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
    dry_array = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
    show_coords = [[] for p in show_particles]
    for i in range(particles.n_particles):
        if x_axis.valid_bin(x_bin[i]) and y_axis.valid_bin(y_bin[i]):
            if water_frac[i] > 0.0:
                if water_frac_array[x_bin[i], y_bin[i]] == 0.0:
                    water_frac_array[x_bin[i], y_bin[i]] = water_frac[i]
                else:
                    water_frac_array[x_bin[i], y_bin[i]] \
                         = min(water_frac_array[x_bin[i], y_bin[i]], water_frac[i])
            else:
                dry_array[x_bin[i], y_bin[i]] = 1.0
            for j in range(len(show_particles)):
                if particles.id[i] == show_particles[j]["id"]:
                    show_coords[j] = [diameter[i], comp_frac[i]]

    if print_info:
        print "%g hours, %s LST, water = %g%% to %g%%" \
              % (env_state.elapsed_time / 3600,
                 time_of_day_string(env_state.elapsed_time
                                    + env_state.start_time_of_day),
                 water_frac_array[water_frac_array > 0.0].min(),
                 water_frac_array.max())

    value = water_frac_array
    if value_max != None:
        value = value / value_max
    value = value.clip(0.0, 1.0)
    value_mask = (value != 0.0)
    value = value * (1 - min_palette_index) + min_palette_index

    stroke_width_array = dry_array * 0.0001

    wet_rects = pmc_histogram_2d_multi([value], x_axis, y_axis,
                                       mask = value_mask)
    dry_rects = pmc_histogram_2d_multi([ones_like(dry_array)],
                                       x_axis, y_axis, mask = dry_array)
    white_rects = pmc_histogram_2d_multi([zeros_like(dry_array)],
                                         x_axis, y_axis, inv_mask = dry_array)
    return (wet_rects, dry_rects, dry_array, show_coords, env_state,
            white_rects)

time_filename_list = get_time_filename_list(netcdf_dir_nc, netcdf_pattern_nc)
count = 0
start_time = module_time.clock()
for [time, filename, output_key] in time_filename_list:
    g = make_1x1_graph_grid(y_axis_label)

    (wet_rects, dry_rects, dry_mask_array, show_coords, env_state,
     white_rects) = get_plot_data(filename, max_val,
                                  print_info = False)

    if len(dry_rects) > 0:
        g.plot(graph.data.points(dry_rects,
                                 xmin = 1, xmax = 2, ymin = 3, ymax = 4,
                                 color = 5),
               styles = [hsb_rect(gray_palette)])

    if len(wet_rects) > 0:
        g.plot(graph.data.points(wet_rects,
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
                         density = 0.8,
                         min_palette_index = min_palette_index,
                         xpos = g.xpos + g.width + grid_h_space,
                         ybottom = g.ypos,
                         ytop = g.ypos + g.height,
                         title = r"water mass frac. $w_{\rm H_2O,all}$ ($\%$)",
                         palette = rainbow_palette,
                         extra_box_color = color.gray(0),
                         extra_box_pattern = False,
                         extra_box_label = "dry")

    out_filename = "%s_%s.pdf" % (out_prefix, output_key)
    g.writePDFfile(out_filename)

    count += 1
    total_count = len(time_filename_list)
    elapsed_time = module_time.clock() - start_time
    total_time = float(total_count) / float(count) * elapsed_time
    remaining_time = total_time - elapsed_time
    print "water, no coag: %d of %d - %.1f s elapsed, %.1f s remaining" \
          % (count, total_count, elapsed_time, remaining_time)
