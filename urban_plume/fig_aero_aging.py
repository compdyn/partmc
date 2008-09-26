#!/usr/bin/env python
# Copyright (C) 2007, 2008 Matthew West
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

out_prefix = "figs/aero_aging"

time_filename_list_wc = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max([time for [time, filename, key] in time_filename_list_wc]) / 60

g = graph.graphxy(
    width = 6.7,
    x = graph.axis.linear(min = 0.,
                          max = max_time_min,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time
                                               = start_time_of_day_min),
                          title = "local standard time (hours:minutes)",
                          painter = grid_painter),
    y = graph.axis.linear(title = r"number fraction"))

plot_data = []
for [time, filename, key] in time_filename_list_wc:
    print filename
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()
    num_den = 1.0 / particles.comp_vol
    critical_ss = particles.kappa_rh(env_state) - 1.0
    total_num_den = num_den.sum()
    aged_num_den = 0.0
    for i in range(num_den.size):
        if critical_ss[i] < 0.05:
            aged_num_den += num_den[i]
    value = aged_num_den / total_num_den
    plot_data.append([time / 60.0, value])

g.plot(
    graph.data.points(plot_data, x = 1, y = 2),
    styles = [graph.style.line()])

out_filename = "%s_color.pdf" % out_prefix
g.writePDFfile(out_filename)

