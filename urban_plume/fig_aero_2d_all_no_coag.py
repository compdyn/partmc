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

y_axis_label = r"$f_{{\rm BC},{\rm all}}$ ($1$)"
out_filename = "figs/aero_2d_all_no_coag.pdf"

def get_plot_data(filename, value_max = None):
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()

    diameter = particles.dry_diameter() * 1e6
    comp_frac = particles.mass(include = ["BC"]) \
                / particles.mass(exclude = ["H2O"]) * 100

    x_axis = pmc_log_axis(min = diameter_axis_min, max = diameter_axis_max,
                          n_bin = 70)
    y_axis = pmc_linear_axis(min = 0, max = 100, n_bin = 100)
    x_bin = x_axis.find(diameter)
    y_bin = y_axis.find(comp_frac)

    num_den_array = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
    for i in range(particles.n_particles):
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

graphs = make_2x2_graph_grid(y_axis_label)
time_filename_list = get_time_filename_list(netcdf_dir_nc, netcdf_pattern_nc)
for (graph_name, time_hour) in times_hour.iteritems():
    time = time_hour * 3600.0
    filename = file_filename_at_time(time_filename_list, time)
    (rects, env_state) = get_plot_data(filename, max_val)
    g = graphs[graph_name]
    g.plot(graph.data.points(rects,
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

c = graphs["c"]
add_canvas_color_bar(c,
                     min = 0.0,
                     max = max_val,
                     title = r"normalized number density $\hat{n}_{\rm BC,dry}(f,D)$ (1)",
                     palette = gray_palette)

c.writePDFfile(out_filename)
print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
