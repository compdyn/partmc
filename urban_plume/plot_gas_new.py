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
netcdf_re = re.compile(r"urban_plume_state_0001_([0-9]{8})\.nc")

gas_history = read_history(gas_state_t, netcdf_dir, netcdf_pattern)

g = graph.graphxy(
    width = 10,
    x = graph.axis.linear(min = 0.,
                          max = 1440,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time = 6 * 60),
                          title = "time (LST)",
			  painter = grid_painter),
    y = graph.axis.log(min = 0.01,
                       max = 200,
                       title = "gas concentration (ppb)",
                       painter = grid_painter),
    key = graph.key.key(pos = "br"))

for i in range(len(gas_species)):
    data_slice = module_copy.deepcopy(data)
    data_slice.reduce([select("gas_species", gas_species[i])])
    #data_slice.scale_dim("time", 1.0/3600)
    g.plot(graph.data.points(data_slice.data_center_list(strip_zero = True),
                             x = 1, y = 2,
                             title = tex_species(gas_species[i])),
	   styles = [graph.style.line(lineattrs = [color_list[i], style.linewidth.THick])])

g.writePDFfile("out/gas.pdf")
