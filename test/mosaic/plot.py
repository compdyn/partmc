#!/usr/bin/env python
# Copyright (C) 2007 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../../tool")
from pmc_data_nc import *
from pmc_pyx import *

aero_species = ["SO4_a", "NO3_a", "NH4_a"]
gas_species = ["NO2", "NH3"]

aero_data = pmc_var(NetCDFFile("out/mosaic_state_0001.nc"),
		    "aero",
		    [sum("radius")])
gas_data = pmc_var(NetCDFFile("out/mosaic_state_0001.nc"),
		   "gas",
		   [])

i_time = aero_data.find_dim_by_name("time")
max_time = max(aero_data.dims[i_time].grid_centers) / 3600.0
g = graph.graphxy(width = 10,
		  x = graph.axis.linear(max = max_time,
					title = "time (hour)",
					painter = grid_painter),
		  y = graph.axis.linear(title = "aerosol volume density (m$^3$/m$^3$)"),
		  y2 = graph.axis.linear(title = "gas concentration (ppb)"),
		  key = graph.key.key(pos = "tl"))

for i in range(len(aero_species)):
    data_slice = module_copy.deepcopy(aero_data)
    data_slice.reduce([select("unit", "vol_den"),
		       select("aero_species", aero_species[i])])
    data_slice.scale_dim("time", 1.0/3600)
    g.plot(graph.data.list(data_slice.data_center_list(),
			   x = 1, y = 2,
			   title = "aerosol %s" % tex_species(aero_species[i])),
	   styles = [graph.style.line(lineattrs = [color_list[i]])])

for i in range(len(gas_species)):
    data_slice = module_copy.deepcopy(gas_data)
    data_slice.reduce([select("gas_species", gas_species[i])])
    data_slice.scale_dim("time", 1.0/3600)
    g.plot(graph.data.list(data_slice.data_center_list(),
			   x = 1, y2 = 2,
			   title = "gas %s" % tex_species(gas_species[i])),
	   styles = [graph.style.line(lineattrs = [color_list[i],
						   style.linestyle.dashed])])

g.writePDFfile("out/mosaic.pdf")

