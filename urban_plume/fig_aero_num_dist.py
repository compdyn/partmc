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

disp_lines = [
    {"time_hour": 1, "coag": True,
     "line_style": 2, "line_thickness": style.linewidth.Thick,
     "line_color": 0, "line_color_style": style.linestyle.solid},
    {"time_hour": 6, "coag": True,
     "line_style": 1, "line_thickness": style.linewidth.Thick,
     "line_color": 1, "line_color_style": style.linestyle.solid},
    {"time_hour": 24, "coag": True,
     "line_style": 0, "line_thickness": style.linewidth.Thick,
     "line_color": 2, "line_color_style": style.linestyle.solid},
    {"time_hour": 24, "coag": False,
     "line_style": 0, "line_thickness": style.linewidth.THIck,
     "line_color": 2, "line_color_style": style.linestyle.dashed},
    ]

print_diams = [0.03, 0.05, 0.07, 0.10]
base_vals = []
new_vals = []
eval_change_time = 24

out_prefix = "figs/aero_num_dist"

x_axis = pmc_log_axis(min = diameter_axis_min, max = diameter_axis_max,
                      n_bin = num_diameter_bins)

time_filename_list_wc = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
time_filename_list_nc = get_time_filename_list(netcdf_dir_nc, netcdf_pattern_nc)

for use_color in [True, False]:
    g = graph.graphxy(
        width = 6.8,
        x = graph.axis.log(min = x_axis.min,
                           max = x_axis.max,
                           title = r"dry diameter $D$ ($\rm \mu m$)",
                           painter = major_grid_painter),
        y = graph.axis.log(min = 1e7,
                           max = 1e11,
                           title = r"number concentration ($\rm m^{-3}$)",
                           painter = major_grid_painter),
        key = graph.key.key(pos = None, hpos = 0.6, vpos = 0))

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

        x_bin = x_axis.find(diameter)

        num_den_array = numpy.zeros([x_axis.n_bin])
        for i in range(particles.n_particles):
            scale = particles.comp_vol[i] * x_axis.grid_size(x_bin[i])
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
        if not use_color:
            for d in print_diams:
                x_bin = x_axis.find([d])[0]
                print "time = %g hours, coag = %s, n(%g) = %g m^{-3}" \
                      % (disp_lines[t]["time_hour"],
                         str(disp_lines[t]["coag"]), d, num_den_array[x_bin])
                if disp_lines[t]["time_hour"] == eval_change_time:
                    if disp_lines[t]["coag"] == False:
                        base_vals.append(num_den_array[x_bin])
                    else:
                        new_vals.append(num_den_array[x_bin])

    if not use_color:
        if len(base_vals) != len(print_diams):
            print "ERROR: something wrong with base_vals"
        if len(new_vals) != len(print_diams):
            print "ERROR: something wrong with new_vals"
        for i in range(len(print_diams)):
            print ("decrease in n(%g) from no-coag to with-coag "
                   "at %g hours = %g%%") \
                   % (print_diams[i], eval_change_time,
                      (base_vals[i] - new_vals[i]) / base_vals[i] * 100)

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    g.writePDFfile(out_filename)
    if not color:
        print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
