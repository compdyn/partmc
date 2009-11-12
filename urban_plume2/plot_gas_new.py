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

gas_species = ["O3", "NO2", "HCHO", "HNO3", "SO2", "NH3"]

netcdf_dir = "out"
netcdf_pattern = r"urban_plume_state_0001_([0-9]{8})\.nc"

gas_state_history = read_history(gas_state_t, netcdf_dir, netcdf_pattern)
env_state = read_any(env_state_t, netcdf_dir, netcdf_pattern)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max([time for [time, gas_state] in gas_state_history]) / 60

g = graph.graphxy(
    width = 10,
    x = graph.axis.linear(min = 0.0,
                          max = max_time_min,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time
                                               = start_time_of_day_min),
                          title = "time (LST)",
			  painter = grid_painter),
    y = graph.axis.log(min = 0.01,
                       max = 200,
                       title = "gas concentration (ppb)",
                       painter = grid_painter),
    key = graph.key.key(pos = "br"))

for i in range(len(gas_species)):
    plot_data = []
    for [time, gas_state] in gas_state_history:
        conc = gas_state.concentration_by_species(gas_species[i])
        if conc > 0.0:
            plot_data.append([time / 60, conc])
    if plot_data != []:
        g.plot(graph.data.points(plot_data, x = 1, y = 2,
                                 title = tex_species(gas_species[i])),
               styles = [graph.style.line(lineattrs = [color_list[i],
                                                       style.linewidth.THick])])
    else:
        print "warning: only zeros for species %s" % gas_species[i]

g.writePDFfile("out/gas.pdf")
