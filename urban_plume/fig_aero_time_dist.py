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

aero_species = ["NO3", "NH4", "OC", "SO4","BC"]
#aero_species = ["BC" ]
#aero_species = ["ARO1", "ARO2", "ALK1", "OLE1" ]
#aero_species = ["API1", "API2", "LIM1", "LIM2" ]
#aero_species = ["OIN" ]

subdir = "withcoag_dry"
if len(sys.argv) > 1:
    subdir = sys.argv[1]

data = pmc_var(NetCDFFile("out/%s/urban_plume_0001.nc" % subdir),
	       "aero",
	       [sum("dry_radius"),
		select("unit", "mass_den")])
data.write_summary(sys.stdout)

data.scale_dim("time", 1.0/60)
data.scale(1e9)

g = graph.graphxy(
    width = 10,
    x = graph.axis.linear(min = 0,
                          max = 1440,
                          title = "local standard time",
#                          parter = graph.axis.parter.linear(tickdists = [6, 3]),
#                          max = max(data.dim_by_name("time").grid_centers),
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time = 6 * 60),
			  painter = grid_painter),
    y = graph.axis.linear(title = "mass density ($\mu$g/m$^3$)",
			  painter = grid_painter),
    key = graph.key.key(pos = "tl"))
print max(data.dim_by_name("time").grid_centers)
for i in range(len(aero_species)):
    data_slice = module_copy.deepcopy(data)
    data_slice.reduce([select("aero_species", aero_species[i])])
    g.plot(graph.data.list(data_slice.data_center_list(),
			   x = 1, y = 2,
			   title = tex_species(aero_species[i])),
	   styles = [graph.style.line(lineattrs = [color_list[i],style.linewidth.THick])])

#data.reduce([sum("aero_species")])
#g.plot(graph.data.list(data.data_center_list(),
 #                      x = 1, y = 2,
 #                      title = "total"),
 #      styles = [graph.style.line(lineattrs
 #                                 = [color_list[len(aero_species)]])])

g.writePDFfile("out/%s/aero_hist.pdf" % subdir)
