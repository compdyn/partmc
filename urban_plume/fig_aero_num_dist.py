#!/usr/bin/env python
# Copyright (C) 2007, 2008 Matthew West
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

times_hour = [1, 6, 24, 24]
with_coag = [True, True, True, False]
line_style_order = [2, 1, 0, 0]
print_diams = [0.03, 0.05, 0.07, 0.10]

out_filename = "figs/aero_num_dist.pdf"

x_axis = pmc_log_axis(min = diameter_axis_min, max = diameter_axis_max,
                      n_bin = 70)

g = graph.graphxy(
    width = 6.8,
    x = graph.axis.log(min = x_axis.min,
		       max = x_axis.max,
		       title = r"dry diameter ($\rm \mu m$)",
		       painter = major_grid_painter),
    y = graph.axis.log(min = 1e7,
                       max = 1e11,
                       title = r"number density ($\rm m^{-3}$)",
		       painter = major_grid_painter),
    key = graph.key.key(pos = None, hpos = 0.6, vpos = 0))

time_filename_list_wc = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
time_filename_list_nc = get_time_filename_list(netcdf_dir_nc, netcdf_pattern_nc)
for t in range(len(times_hour)):
    if with_coag[t]:
        filename = file_filename_at_time(time_filename_list_wc,
                                         times_hour[t] * 3600)
    else:
        filename = file_filename_at_time(time_filename_list_nc,
                                         times_hour[t] * 3600)
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    ncf.close()

    diameter = particles.dry_diameter() * 1e6

    x_bin = x_axis.find(diameter)

    num_den_array = numpy.zeros([x_axis.n_bin])
    for i in range(particles.n_particles):
        scale = particles.comp_vol[i] * x_axis.grid_size(x_bin[i])
        num_den_array[x_bin[i]] += 1.0 / scale

    plot_data = [[x_axis.center(i), num_den_array[i]]
                 for i in range(x_axis.n_bin) if num_den_array[i] > 0.0]

    if times_hour[t] == 1:
        title = "1 hour"
    else:
        title = "%g hours" % times_hour[t]
    if not with_coag[t]:
        title = "%s, no coag" % title
    if t < 3:
        thickness = style.linewidth.Thick
    else:
        thickness = style.linewidth.THIck
    g.plot(graph.data.points(plot_data, x = 1, y = 2, title = title),
           styles = [graph.style.line(lineattrs
                                      = [line_style_list[line_style_order[t]],
                                         thickness])])
    for d in print_diams:
        x_bin = x_axis.find([d])[0]
        print "time = %g hours, coag = %s, n(%g) = %g m^{-3}" \
              % (times_hour[t], str(with_coag[t]), d, num_den_array[x_bin])

g.writePDFfile(out_filename)
print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
