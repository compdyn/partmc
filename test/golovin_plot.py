#!/usr/bin/env python
# Copyright (C) 2007 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
sys.path.append(os.path.expanduser("~/.python"))
from pyx import *
from Scientific.IO.NetCDF import *

times_min = [0, 5, 10]
times_sec = [t * 60 for t in times_min]

mc_data = pmc_var(NetCDFFile("out/golovin_mc_state_0001.nc"),
		  "aero",
		  [sum("aero_species")])
exact_data = pmc_var(NetCDFFile("out/golovin_exact_0001.nc"),
		     "aero",
		     [sum("aero_species")])

g_num = graph.graphxy(width = 10,
		      x = graph.axis.log(min = 1e-7,
					 max = 1e-3,
					 title = "radius (m)",
					 painter = grid_painter),
		      y = graph.axis.log(min = 1e5,
					 max = 1e10,
					 title = "number density (\#/m$^3$)",
					 painter = grid_painter),
		      key = graph.key.key(pos = "tr"))
g_vol = graph.graphxy(width = 10,
		      x = graph.axis.log(min = 1e-7,
					 max = 1e-3,
					 title = "radius (m)",
					 painter = grid_painter),
		      y = graph.axis.log(min = 1e-13,
					 max = 1e-5,
					 title = "volume density (m$^3$/m$^3$)",
					 painter = grid_painter),
		      key = graph.key.key(pos = "br"))

for i in range(len(times_sec)):
    data_slice = module_copy.deepcopy(mc_data)
    data_slice.reduce([select("unit", "num_den"),
		       select("time", times_sec[i])])
    g_num.plot(graph.data.list(data_slice.data_center_list(strip_zero = True),
			       x = 1, y = 2,
			       title = "%g mins MC" % times_min[i]),
	       styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
					    size = 0.05,
					    symbolattrs = [color_list[i]])])

    data_slice = module_copy.deepcopy(mc_data)
    data_slice.reduce([select("unit", "vol_den"),
		       select("time", times_sec[i])])
    g_vol.plot(graph.data.list(data_slice.data_center_list(strip_zero = True),
			       x = 1, y = 2,
			       title = "%g mins MC" % times_min[i]),
	       styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
					    size = 0.05,
					    symbolattrs = [color_list[i]])])

    data_slice = module_copy.deepcopy(exact_data)
    data_slice.reduce([select("unit", "num_den"),
		       select("time", times_sec[i])])
    g_num.plot(graph.data.list(data_slice.data_center_list(strip_zero = True),
			       x = 1, y = 2,
			       title = "%g mins exact" % times_min[i]),
	       styles = [graph.style.line(lineattrs = [color_list[i]])])
    
    data_slice = module_copy.deepcopy(exact_data)
    data_slice.reduce([select("unit", "vol_den"),
		       select("time", times_sec[i])])
    g_vol.plot(graph.data.list(data_slice.data_center_list(strip_zero = True),
			       x = 1, y = 2,
			       title = "%g mins exact" % times_min[i]),
	       styles = [graph.style.line(lineattrs = [color_list[i]])])
    
g_num.writePDFfile("out/golovin_num.pdf")
g_vol.writePDFfile("out/golovin_vol.pdf")
