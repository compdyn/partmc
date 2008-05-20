#!/usr/bin/env python
# Copyright (C) 2007-2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../../tool")
from pmc_data_nc import *
from pmc_pyx import *

times_min = [0, 5, 10]
times_sec = [t * 60 for t in times_min]

mc_data = pmc_var(NetCDFFile("out/golovin_mc_avg.nc"),
		  "aero",
		  [sum("aero_species")])
exact_data = pmc_var(NetCDFFile("out/golovin_exact_0001.nc"),
		     "aero",
		     [sum("aero_species")])

mc_data.scale_dim("radius", 2e6)
exact_data.scale_dim("radius",2e6)

g_num = graph.graphxy(width = 6.6,
		      x = graph.axis.log(min = 1e-1,
					 max = 1e3,
					 title = r"diameter ($\rm \mu m$)",
					 painter = major_grid_painter),
		      y = graph.axis.log(min = 1e5,
					 max = 1e10,
					 title = r"number density ($\rm m^{-3}$)",
					 painter = major_grid_painter))
g_vol = graph.graphxy(width = 6.6,
		      x = graph.axis.log(min = 1e-1,
					 max = 1e3,
					 title = "diameter ($\mu$m)",
					 painter = grid_painter),
		      y = graph.axis.log(min = 1e-13,
					 max = 1e-5,
					 title = "volume density (m$^3$/m$^3$)",
					 painter = grid_painter))

for i in range(len(times_sec)):
    data_slice = module_copy.deepcopy(mc_data)
    data_slice.reduce([select("unit", "num_den"),
		       select("time", times_sec[i])])
    g_num.plot(graph.data.points(data_slice.data_center_list(strip_zero = True),
			       x = 1, y = 2,
			       title = "%g mins MC" % times_min[i]),
	       styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
					    size = 0.05,
					    symbolattrs = [color.grey.black])])

    data_slice = module_copy.deepcopy(exact_data)
    data_slice.reduce([select("unit", "num_den"),
		       select("time", times_sec[i])])
    g_num.plot(graph.data.points(data_slice.data_center_list(strip_zero = True),
			       x = 1, y = 2,
			       title = "%g mins exact" % times_min[i]),
	       styles = [graph.style.line(lineattrs = [color.grey.black])])

(x, y) = g_num.pos(37, 4e8)
g_num.text(x, y, "0 mins", [text.halign.left, text.valign.bottom])
(x, y) = g_num.pos(80, 3e7)
g_num.text(x, y, "5 mins", [text.halign.left, text.valign.bottom])
(x, y) = g_num.pos(160, 3e6)
g_num.text(x, y, "10 mins", [text.halign.left, text.valign.bottom])

g_num.writePDFfile("out/golovin.pdf")

print "figure height = %.1f cm" % unit.tocm(g_num.bbox().height())
print "figure width = %.1f cm" % unit.tocm(g_num.bbox().width())
