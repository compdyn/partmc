#!/usr/bin/env python
# Copyright (C) 2007-2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *

from fig_helper import *

y_axis_label = r"BC fraction $f_{{\rm BC},{\rm dry}}$ ($1$)"
out_filename = "figs/aero_2d_water.pdf"

max_val = 50
min_palette_index = 0

fill_pattern = hash_pattern(n_lines = 7,
                            line_attributes = [style.linewidth.normal])

def get_plot_data(filename, value_max = None):
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
                          n_bin = 70)
    y_axis = pmc_linear_axis(min = 0, max = 100, n_bin = 100)
    x_bin = x_axis.find(diameter)
    y_bin = y_axis.find(comp_frac)

    water_frac_array = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
    dry_array = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
    show_coords = [[] for p in show_particles]
    for i in range(particles.n_particles):
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
    return (wet_rects, dry_rects, dry_array, show_coords, env_state)

graphs = make_2x2_graph_grid(y_axis_label)
time_filename_list = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
for (graph_name, time_hour) in times_hour.iteritems():
    time = time_hour * 3600.0
    filename = file_filename_at_time(time_filename_list, time)
    (wet_rects, dry_rects, dry_mask_array, show_coords, env_state) \
                = get_plot_data(filename, max_val)
    g = graphs[graph_name]
    if len(dry_rects) > 0:
        g.plot(graph.data.points(dry_rects,
                                 xmin = 1, xmax = 2, ymin = 3, ymax = 4,
                                 color = 5),
               styles = [hsb_rect(gray_palette,
                                  fill_pattern = fill_pattern,
                                  do_stroke = False)])
    g.plot(graph.data.points(wet_rects,
                             xmin = 1, xmax = 2, ymin = 3, ymax = 4,
                             color = 5),
           styles = [hsb_rect(gray_palette)])

    g.dolayout()
    for axisname in ["x", "y"]:
        for t in g.axes[axisname].data.ticks:
            if t.ticklevel is not None:
                g.stroke(g.axes[axisname].positioner.vgridpath(t.temp_v),
                         [style.linestyle.dotted])
    g.dodata()
    g.doaxes()

    write_time(g, env_state)
    for i in range(len(show_particles)):
        if len(show_coords[i]) > 0:
            label_point(g, show_coords[i][0], show_coords[i][1],
                        show_particles[i]["label pos"][0],
                        show_particles[i]["label pos"][1],
                        show_particles[i]["label"])

c = graphs["c"]
add_canvas_color_bar(c,
                     min = 0.0,
                     max = max_val,
                     min_palette_index = min_palette_index,
                     title = r"water fraction $f_{\rm H_2O,all}$ (1)",
                     texter = graph.axis.texter.decimal(suffix = r"\%"),
                     palette = gray_palette,
                     extra_box_value = 0.0,
                     extra_box_pattern = fill_pattern,
                     extra_box_label = "dry")

c.writePDFfile(out_filename)
print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
