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

aero_species = [["", ["NO3"]],
                ["", ["NH4"]],
                ["", ["OC"]],
                ["", ["SO4"]],
                ["", ["BC"]],
                ["SOA", ["ARO1", "ARO2", "ALK1", "OLE1"]],
                ["", ["H2O"]],
                ]

data = pmc_var(NetCDFFile("out/urban_plume_0.5_3am_0001.nc"),
	       "aero",
	       [sum("dry_radius"),
		select("unit", "mass_den")])

data.scale_dim("time", 1.0/60)
data.scale(1e9)

g = graph.graphxy(
    width = 10,
    x = graph.axis.linear(min = 0.,
                          max = 1440,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time = 3 * 60),
                          title = "time (LST)",
			  painter = grid_painter),
    y = graph.axis.log(title = r"mass density ($\rm \mu g \, m^{-3}$)",
                       painter = grid_painter),
    key = graph.key.key(pos = "br"))

for i in range(len(aero_species)):
    data_slice = module_copy.deepcopy(data)
    data_slice.reduce([sum("aero_species", only = aero_species[i][1])])
    if aero_species[i][0] != "":
        title = aero_species[i][0]
    else:
        title = tex_species(aero_species[i][1][0])
    g.plot(graph.data.points(data_slice.data_center_list(strip_zero = True),
                             x = 1, y = 2,
                             title = title),
	   styles = [graph.style.line(lineattrs = [color_list[i], style.linewidth.Thick])])

g.writePDFfile("out/aero_bulk.pdf")
