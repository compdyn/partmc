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

aero_species = ["SO4", "NO3", "NH4", "OC"]
#aero_species = ["BC" ]
#aero_species = ["ARO1", "ARO2", "ALK1", "OLE1" ]
#aero_species = ["API1", "API2", "LIM1", "LIM2" ]
#aero_species = ["OIN" ]

data = pmc_var(NetCDFFile("out/testcase_withcoag/urban_plume_state_0001.nc"),
	       "aero",
	       [sum("radius"),
		select("unit", "mass_den")])
data.write_summary(sys.stdout)

data.scale_dim("time", 1.0/3600)
data.scale(1e9)

g = graph.graphxy(
    width = 10,
    x = graph.axis.linear(title = "time (hour)",
			  painter = grid_painter),
    y = graph.axis.linear(title = "mass density ($\mu$g/m$^3$)",
			  painter = grid_painter),
    key = graph.key.key(pos = "tr"))

for i in range(len(aero_species)):
    data_slice = module_copy.deepcopy(data)
    data_slice.reduce([select("aero_species", aero_species[i])])
    print data_slice.data
    g.plot(graph.data.list(data_slice.data_center_list(),
			   x = 1, y = 2,
			   title = tex_species(aero_species[i])),
	   styles = [graph.style.line(lineattrs = [color_list[i]])])

g.writePDFfile("out/testcase_withcoag/aero_hist_check.pdf")
