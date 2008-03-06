#!/usr/bin/env python
# Copyright (C) 2007 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *

#env = ["NO2",  "O3", "HNO3", "SO2", "NH3", "NO"]

data = pmc_var(NetCDFFile("out/testcase_nocoag/urban_plume_state_0001.nc"),
	       "env",
	       [])

data.write_summary(sys.stdout)

data.reduce([select("env", "rel_humid")])
print data.data
#g = graph.graphxy(
#    width = 10,
#    x = graph.axis.linear(title = "time (hour)",
#			  painter = grid_painter),
#    y = graph.axis.linear(title = "gas concentration (ppb)",
#			  painter = grid_painter),
#    key = graph.key.key(pos = "tr"))

#for i in range(len(gas_species)):
#    data_slice = module_copy.deepcopy(data)
#    data_slice.reduce([select("gas_species", gas_species[i])])
#    data_slice.scale_dim("time", 1.0/3600)
#    g.plot(graph.data.list(data_slice.data_center_list(),
#			   x = 1, y = 2,
#			   title = tex_species(gas_species[i])),
#	   styles = [graph.style.line(lineattrs = [color_list[i]])])

#g.writePDFfile("out/testcase_nocoag/gas_dilut_anorg_4.pdf")
