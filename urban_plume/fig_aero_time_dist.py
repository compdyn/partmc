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

aero_species = [["NO3"],
                ["NH4"],
                ["OC"],
                ["SO4"],
                ["BC"],
                ["ARO1", "ARO2", "ALK1", "OLE1" ]]
line_style_order = [4, 5, 0, 1, 2, 3]
key_names = [ None, None, None, None, None, "SOA"]  # None means use default

data = pmc_var(NetCDFFile("out/urban_plume_with_coag_0001.nc"),
	       "aero",
	       [sum("dry_radius"),
		select("unit", "mass_den")])
#data.write_summary(sys.stdout)

data.scale_dim("time", 1.0/60)
data.scale(1e9)

g = graph.graphxy(
    width = 6.7,
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
    key = graph.key.key(pos = "tc",
                        vinside = 0,
                        columns = 3,
                        symbolwidth = 0.75 * unit.v_cm,
                        keyattrs = [deco.stroked, deco.filled([color.gray.white])]))
# print max(data.dim_by_name("time").grid_centers)
for i in range(len(aero_species)):
    data_slice = module_copy.deepcopy(data)
    data_slice.reduce([sum("aero_species", only = aero_species[i])])
    if key_names[i] == None:
        title_name =  tex_species(aero_species[i][0])
    else:
        title_name = key_names[i]
        
    g.plot(graph.data.points(data_slice.data_center_list(),
			   x = 1, y = 2,
			   title = title_name),
	   styles = [graph.style.line(lineattrs = [line_style_list[line_style_order[i]],style.linewidth.Thick])])

#data.reduce([sum("aero_species")])
#g.plot(graph.data.points(data.data_center_list(),
 #                      x = 1, y = 2,
 #                      title = "total"),
 #      styles = [graph.style.line(lineattrs
 #                                 = [color_list[len(aero_species)]])])

g.writePDFfile("figs/aero_time_dist.pdf")
print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
