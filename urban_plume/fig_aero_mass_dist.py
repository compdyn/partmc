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

out_filename = "figs/aero_mass_dist.pdf"

x_axis = pmc_log_axis(min = diameter_axis_min, max = diameter_axis_max,
                      n_bin = 70)

c = canvas.canvas()

g2 = c.insert(graph.graphxy(
    width = 6.8,
    x = graph.axis.log(min = x_axis.min,
		       max = x_axis.max,
		       title = r"dry diameter ($\rm \mu m$)",
		       painter = major_grid_painter),
    y = graph.axis.log(min = 1e-6,
                       max = 1e2,
                       title = r"mass density ($\rm \mu g \, m^{-3}$)",
		       painter = major_grid_painter),
    key = graph.key.key(pos = None, hpos = 0.8, vpos = 0)))

g1 = c.insert(graph.graphxy(
    width = 6.8,
    ypos = g2.height + 0.5,
    x = graph.axis.linkedaxis(g2.axes["x"],
                              painter = graph.axis.painter.linked(gridattrs = [attr.changelist([style.linestyle.dotted, None])])),
    y = graph.axis.log(min = 1e-6,
                       max = 1e2,
                       title = r"mass density ($\rm \mu g \, m^{-3}$)",
		       painter = major_grid_painter),
    key = graph.key.key(pos = "br")))

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
    so4_mass = particles.mass(include = ["SO4"]) * 1e9
    bc_mass = particles.mass(include = ["BC"]) * 1e9

    x_axis = pmc_log_axis(min = 1e-2, max = 2, n_bin = 70)
    x_bin = x_axis.find(diameter)

    so4_mass_array = numpy.zeros([x_axis.n_bin])
    bc_mass_array = numpy.zeros([x_axis.n_bin])
    for i in range(particles.n_particles):
        scale = particles.comp_vol[i] * x_axis.grid_size(x_bin[i])
        so4_mass_array[x_bin[i]] += so4_mass[i] / scale
        bc_mass_array[x_bin[i]] += bc_mass[i] / scale

    so4_plot_data = [[x_axis.center(i), so4_mass_array[i]]
                     for i in range(x_axis.n_bin) if so4_mass_array[i] > 0.0]
    bc_plot_data = [[x_axis.center(i), bc_mass_array[i]]
                    for i in range(x_axis.n_bin) if so4_mass_array[i] > 0.0]

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
    g1.plot(graph.data.points(so4_plot_data, x = 1, y = 2, title = title),
            styles = [graph.style.line(lineattrs
                                       = [line_style_list[line_style_order[t]],
                                          thickness])])
    g2.plot(graph.data.points(bc_plot_data, x = 1, y = 2, title = title),
            styles = [graph.style.line(lineattrs
                                       = [line_style_list[line_style_order[t]],
                                          thickness])])

write_text_inside(g1, r"$\rm SO_4$")
write_text_inside(g2, r"$\rm BC$")

c.writePDFfile(out_filename)
print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
