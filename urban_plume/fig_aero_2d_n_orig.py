#!/usr/bin/env python
# Copyright (C) 2007-2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
from fig_helper import *

time_hour = 24
max_n_coags = 30

y_axis_label = r"$f_{{\rm BC},{\rm all}}$ ($1$)"
out_filename = "figs/aero_2d_n_orig.pdf"

netcdf_dir = "out"
netcdf_pattern = r"urban_plume_state_0001_([0-9]{8})\.nc"


g = graph.graphxy(
    width = 6.9,
    x = graph.axis.log(min = 1e-3,
                       max = 1e+0,
                       title = r'dry radius ($\rm \mu m$)',
                       painter = grid_painter),
    y = graph.axis.linear(min = 0,
                          max = max_n_coags,
                          title = 'number of coagulation events',
                          painter = grid_painter))

time_filename_list = get_time_filename_list(netcdf_dir, netcdf_pattern)
time = time_hour * 3600.0
filename = file_filename_at_time(time_filename_list, time)
ncf = NetCDFFile(filename)
particles = aero_particle_array_t(ncf)
ncf.close()

diameter = particles.dry_diameter() * 1e6

x_axis = pmc_log_axis(min = 1e-2, max = 2, n_bin = 70)
y_axis = pmc_linear_axis(min = 0, max = max_n_coags, n_bin = max_n_coags)
x_bin = x_axis.find(diameter)

num_den_array = numpy.zeros([x_axis.n_bin, max_n_coags])
for i in range(particles.n_particles):
    scale = particles.comp_vol[i] / x_axis.grid_size(x_bin[i])
    n_coags = particles.n_orig_part[i] - 1
    n_coags = min(n_coags, max_n_coags - 1)
    num_den_array[x_bin[i], n_coags] += 1.0 / scale

value = num_den_array
value_max = value.max()
if value_max > 0.0:
    value = value / value_max
value = value.clip(0.0, 1.0)

rects = pmc_histogram_2d_multi([value],
                               x_axis, y_axis)
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

suffix = "s"
if time_hour == 1:
    suffix = ""
boxed_text(g, 0.04, 0.9, "%d hour%s" % (time_hour, suffix))
for i in range(len(show_particles)):
    if len(show_coords[i]) > 0:
        label_point(g, show_coords[i][0], show_coords[i][1],
                    show_particles[i][1][0], show_particles[i][1][1],
                    show_particles[i][2])

add_horiz_color_bar(g,
                    min = 0.0,
                    max = value_max,
                    title = r"number density ($\rm m^{-3}$)",
                    palette = gray_palette,
                    bar_offset = 0.6)

g.writePDFfile(out_filename)
print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
