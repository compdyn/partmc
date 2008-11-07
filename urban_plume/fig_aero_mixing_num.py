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

y_min = 1e6
y_max = 1e11

num_bc_mixing_bins = 4
title_list = [r"0--25\% soot",
              r"25--50\% soot",
              r"50--75\% soot",
              r"75--100\% soot",
              ]
if num_bc_mixing_bins != len(title_list):
    raise Exception("mismatch")

out_prefix = "figs/aero_mixing_num"

y_axis_label = r"number concentration ($\rm m^{-3}$)"

def get_plot_data(filename):
    print filename
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()

    diameter = particles.dry_diameter() * 1e6
    comp_frac = particles.mass(include = ["BC"]) \
                / particles.mass(exclude = ["H2O"]) * 100

    x_axis = pmc_log_axis(min = diameter_axis_min, max = diameter_axis_max,
                          n_bin = num_diameter_bins)
    y_axis = pmc_linear_axis(min = 0, max = 100, n_bin = num_bc_mixing_bins)
    x_bin = x_axis.find(diameter)
    # hack to avoid landing just around the integer boundaries
    comp_frac *= (1.0 + 1e-12)
    y_bin = y_axis.find(comp_frac)

    data_list = [numpy.zeros([x_axis.n_bin])
                 for i in range(num_bc_mixing_bins)]
    for i in range(particles.n_particles):
        scale = particles.comp_vol[i] * x_axis.grid_size(x_bin[i])
        data_list[y_bin[i]][x_bin[i]] += 1.0 / scale

    plot_data_list = [[[x_axis.center(j), data_list[i][j]]
                       for j in range(x_axis.n_bin)
                       if data_list[i][j] > 0.0]
                      for i in range(num_bc_mixing_bins)]
    return (plot_data_list, env_state)

print "***********************************************************************"
print "*                                                                     *"
print "*       WARNING    !!!!!!!    WARNING    !!!!!!!    WARNING           *"
print "*                                                                     *"
print "*                                                                     *"
print "*                      Should be using WC, not NC!                    *"
print "*                                                                     *"
print "*                                                                     *"
print "*       WARNING    !!!!!!!    WARNING    !!!!!!!    WARNING           *"
print "*                                                                     *"
print "***********************************************************************"
time_filename_list = get_time_filename_list(netcdf_dir_nc, netcdf_pattern_nc)
for use_color in [True, False]:
    graphs = make_2x2_graph_grid(y_axis_label, y_min = y_min, y_max = y_max,
                                 with_percent = False, y_log = True,
                                 with_key = True)
    for (graph_name, time_hour) in times_hour.iteritems():
        time = time_hour * 3600.0
        filename = file_filename_at_time(time_filename_list, time)
        (plot_data_list, env_state) = get_plot_data(filename)
        g = graphs[graph_name]
        for i in range(len(plot_data_list)):
            if use_color:
                attrs = [color_list[i], style.linewidth.Thick]
            else:
                attrs = [line_style_list[i], style.linewidth.Thick]
            if len(plot_data_list[i]) > 0:
                g.plot(graph.data.points(plot_data_list[i],
                                         x = 1, y = 2, title = title_list[i]),
                       styles = [graph.style.line(lineattrs = attrs)])

        write_time(g, env_state)

    c = graphs["c"]
    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    c.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
