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
sys.path.append(".")
from fig_helper import *

disp_lines = [
    {"time_hour": 1, "coag": True,
     "line_style": 2, "line_thickness": style.linewidth.Thick,
     "line_color": 0, "line_color_style": style.linestyle.solid},
    {"time_hour": 5, "coag": True,
     "line_style": 1, "line_thickness": style.linewidth.Thick,
     "line_color": 1, "line_color_style": style.linestyle.solid},
    {"time_hour": 7, "coag": True,
     "line_style": 3, "line_thickness": style.linewidth.Thick,
     "line_color": 3, "line_color_style": style.linestyle.solid},
    {"time_hour": 24, "coag": True,
     "line_style": 0, "line_thickness": style.linewidth.Thick,
     "line_color": 2, "line_color_style": style.linestyle.solid},
    {"time_hour": 24, "coag": False,
     "line_style": 0, "line_thickness": style.linewidth.THIck,
     "line_color": 2, "line_color_style": style.linestyle.dashed},
    ]

out_prefix = "figs/aero_mixing_total"

x_axis = pmc_linear_axis(min = bc_axis_min, max = bc_axis_max,
                         n_bin = num_bc_bins)

time_filename_list_wc = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
time_filename_list_nc = get_time_filename_list(netcdf_dir_nc, netcdf_pattern_nc)

for use_color in [True, False]:
    g = graph.graphxy(
        width = 6.8,
        x = graph.axis.linear(min = x_axis.min,
                              max = x_axis.max,
                              title = r"BC dry mass fraction $w_{{\rm BC},{\rm dry}}$ ($1$)",
                              texter = graph.axis.texter.decimal(suffix = r"\%"),
                              painter = grid_painter),
        y = graph.axis.log(min = 1e7,
                           max = 1e12,
                           title = r"number concentration $n_{\rm BC,dry}(w)$ ($\rm m^{-3}$)",
                           painter = major_grid_painter),
        key = graph.key.key(vinside = 0, columns = 2))

    for t in range(len(disp_lines)):
        if disp_lines[t]["coag"]:
            filename = file_filename_at_time(time_filename_list_wc,
                                             disp_lines[t]["time_hour"] * 3600)
        else:
            filename = file_filename_at_time(time_filename_list_nc,
                                             disp_lines[t]["time_hour"] * 3600)
        ncf = NetCDFFile(filename)
        particles = aero_particle_array_t(ncf)
        ncf.close()

        diameter = particles.dry_diameter() * 1e6
        comp_frac = particles.mass(include = ["BC"]) \
            / particles.mass(exclude = ["H2O"]) * 100

        # hack to avoid landing just around the integer boundaries
        comp_frac *= (1.0 + 1e-12)
        x_bin = x_axis.find(comp_frac)

        num_den_array = numpy.zeros([x_axis.n_bin])
        for i in range(particles.n_particles):
            scale = particles.comp_vol[i] * (x_axis.grid_size(x_bin[i]) / 100)
            num_den_array[x_bin[i]] += 1.0 / scale

        plot_data = [[x_axis.center(i), num_den_array[i]]
                     for i in range(x_axis.n_bin) if num_den_array[i] > 0.0]

        if disp_lines[t]["time_hour"] == 1:
            title = "1 hour"
        else:
            title = "%g hours" % disp_lines[t]["time_hour"]
        if not disp_lines[t]["coag"]:
            title = "%s, no coag" % title

        if use_color:
            attrs = [color_list[disp_lines[t]["line_color"]],
                     disp_lines[t]["line_color_style"],
                     style.linewidth.Thick]
        else:
            attrs = [line_style_list[disp_lines[t]["line_style"]],
                     disp_lines[t]["line_thickness"]]
            
        g.plot(graph.data.points(plot_data, x = 1, y = 2, title = title),
               styles = [
            graph.style.line(lineattrs = attrs)])

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    g.writePDFfile(out_filename)
    if not color:
        print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
