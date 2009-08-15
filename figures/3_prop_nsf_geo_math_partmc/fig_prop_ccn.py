#!/usr/bin/env python
# Copyright (C) 2007, 2008, 2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
sys.path.append(".")
from fig_helper import *

out_prefix = "figs_prop/ccn"

critical_ss = 0.6

time_filename_list = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max([time for [time, filename, key] in time_filename_list]) / 60

const = load_constants("../src/constants.f90")

plot_data_resolved = []
plot_data_binned = []
plot_data_ratio = []
for [time, filename, key] in time_filename_list:
    print filename
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()

    diameter_axis = pmc_log_axis(min = diameter_axis_min, max = diameter_axis_max,
                                 n_bin = num_diameter_bins)
    averaged_particles = particles.bin_average(diameter_axis)

    critical_ss_resolved = (particles.kappa_rh(env_state, const) - 1.0) * 100.0
    critical_ss_binned = (averaged_particles.kappa_rh(env_state, const) - 1.0) * 100.0

    bc_mass_resolved = particles.mass(include = ["BC"])
    bc_mass_binned = averaged_particles.mass(include = ["BC"])

    num_den_ccn_resolved = 0.0
    num_den_total_resolved = 0.0
    bc_mass_den_ccn_resolved = 0.0
    bc_mass_den_total_resolved = 0.0
    for i in range(particles.n_particles):
        if bc_mass_resolved[i] > 0:
            num_den_total_resolved += 1.0 / particles.comp_vol[i]
            bc_mass_den_total_resolved += bc_mass_resolved[i] / particles.comp_vol[i]
            if critical_ss_resolved[i] < critical_ss:
                num_den_ccn_resolved += 1.0 / particles.comp_vol[i]
                bc_mass_den_ccn_resolved += bc_mass_resolved[i] / particles.comp_vol[i]
    if num_den_total_resolved == 0:
        num_den_total_resolved = 1
    if bc_mass_den_total_resolved == 0:
        bc_mass_den_total_resolved = 1

    num_den_ccn_binned = 0.0
    num_den_total_binned = 0.0
    bc_mass_den_ccn_binned = 0.0
    bc_mass_den_total_binned = 0.0
    for i in range(averaged_particles.n_particles):
        if bc_mass_binned[i] > 0:
            num_den_total_binned += 1.0 / averaged_particles.comp_vol[i]
            bc_mass_den_total_binned += bc_mass_binned[i] / averaged_particles.comp_vol[i]
            if critical_ss_binned[i] < critical_ss:
                num_den_ccn_binned += 1.0 / averaged_particles.comp_vol[i]
                bc_mass_den_ccn_binned += bc_mass_binned[i] / averaged_particles.comp_vol[i]
    if num_den_total_binned == 0:
        num_den_total_binned = 1
    if bc_mass_den_total_binned == 0:
        bc_mass_den_total_binned = 1

    resolved_value = 100 - num_den_ccn_resolved / num_den_total_resolved * 100
    binned_value = 100 - num_den_ccn_binned / num_den_total_binned * 100
    if binned_value > 0:
        ratio = resolved_value / binned_value
    else:
        ratio = 0
    plot_data_resolved.append([time / 60, resolved_value])
    plot_data_binned.append([time / 60, binned_value])
    plot_data_ratio.append([time / 60, ratio])
    print ratio
    #plot_data_resolved.append([time / 60, bc_mass_den_ccn_resolved / bc_mass_den_total_resolved * 100])
    #plot_data_binned.append([time / 60, bc_mass_den_ccn_binned / bc_mass_den_total_binned * 100])

print max([v for [t,v] in plot_data_ratio])

data_array = zeros([len(plot_data_resolved), 4])
for i in range(len(plot_data_resolved)):
    data_array[i,0] = plot_data_resolved[i][0]
    data_array[i,1] = plot_data_resolved[i][1]
    data_array[i,2] = plot_data_binned[i][1]
    data_array[i,3] = plot_data_ratio[i][1]
savetxt("%s.txt" % out_prefix, data_array)

for use_color in [True, False]:
    g = graph.graphxy(
        width = 6.8,
        x = graph.axis.linear(min = 0.,
                              max = max_time_min,
                              parter = graph.axis.parter.linear(tickdists
                                                                = [6 * 60, 3 * 60]),
                              texter = time_of_day(base_time
                                                   = start_time_of_day_min),
                              title = "local standard time (LST) (hours:minutes)",
                              painter = grid_painter),
        y = graph.axis.linear(min = 0.0,
                              max = 100.0,
                              title = r"unscavenged fraction (\%)",
                              painter = grid_painter))

    g.doaxes()

    if use_color:
        style_attr_resolved = color_list[1]
        style_attr_binned = color_list[2]
    else:
        style_attr_resolved = line_style_list[0]
        style_attr_binned = line_style_list[2]
    g.plot(
        graph.data.points(plot_data_resolved, x = 1, y = 2),
        styles = [graph.style.line(lineattrs
                                   = [style_attr_resolved,
                                      style.linewidth.Thick])])
    g.plot(
        graph.data.points(plot_data_binned, x = 1, y = 2),
        styles = [graph.style.line(lineattrs
                                   = [style_attr_binned,
                                      style.linewidth.Thick])])

    label_plot_line_boxed(g, plot_data_resolved, 4 * 60,
                          "particle-resolved", [0, 0])
    label_plot_line_boxed(g, plot_data_binned, 5 * 60,
                          "sectional", [1, 0])

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    g.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
