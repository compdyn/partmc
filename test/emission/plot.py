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

mc_data = pmc_var(NetCDFFile("out/emission_mc_state_0001.nc"),
		  "aero",
		  [sum("aero_species")])
exact_data = pmc_var(NetCDFFile("out/emission_exact_0001.nc"),
		     "aero",
		     [sum("aero_species")])
sect_data = pmc_var(NetCDFFile("out/emission_sect_0001.nc"),
		    "aero",
		    [sum("aero_species")])

######################################################################

times_min = [0, 1, 10]
times_sec = [t * 60 for t in times_min]

g = graph.graphxy(width = 10,
		  x = graph.axis.log(min = 1e-5,
				     max = 2e-4,
				     title = "radius (m)",
				     painter = grid_painter),
		  y = graph.axis.linear(min = 0,
					max = 1.5e10,
					title = "number density (\#/m$^3$)",
					painter = grid_painter),
		  key = graph.key.key(pos = "tr"))

for i in range(len(times_sec)):
    data_slice = module_copy.deepcopy(mc_data)
    data_slice.reduce([select("unit", "num_den"),
		       select("time", times_sec[i])])
    g.plot(graph.data.list(data_slice.data_center_list(),
			   x = 1, y = 2,
			   title = "%g mins MC" % times_min[i]),
	   styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
					size = 0.15,
					symbolattrs = [color_list[i]])])

    data_slice = module_copy.deepcopy(exact_data)
    data_slice.reduce([select("unit", "num_den"),
		       select("time", times_sec[i])])
    g.plot(graph.data.list(data_slice.data_center_list(),
			   x = 1, y = 2,
			   title = "%g mins exact" % times_min[i]),
	   styles = [graph.style.line(lineattrs = [color_list[i]])])

    data_slice = module_copy.deepcopy(sect_data)
    data_slice.reduce([select("unit", "num_den"),
		       select("time", times_sec[i])])
    g.plot(graph.data.list(data_slice.data_center_list(),
			   x = 1, y = 2,
			   title = "%g mins sect" % times_min[i]),
	   styles = [graph.style.symbol(symbol = graph.style.symbol.cross,
					size = 0.15,
					symbolattrs = [color_list[i]])])
    
g.writePDFfile("out/emission_dist.pdf")

######################################################################

samples = [[1.5e-5, "initial"],
	   [3.0e-5, "emission"],
	   [5.0e-5, "background"]]

g = graph.graphxy(width = 10,
		  x = graph.axis.linear(title = "time (min)",
					painter = grid_painter),
		  y = graph.axis.linear(min = 0,
					max = 1.5e10,
					title = "number density (\#/m$^3$)",
					painter = grid_painter),
		  key = graph.key.key(pos = "br"))

for i in range(len(samples)):
    data_slice = module_copy.deepcopy(mc_data)
    data_slice.reduce([select("unit", "num_den"),
		       select("radius", samples[i][0])])
    data_slice.scale_dim("time", 1.0/60)
    g.plot(graph.data.list(data_slice.data_center_list(),
			   x = 1, y = 2,
			   title = "%s (MC)" % samples[i][1]),
	   styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
					size = 0.05,
					symbolattrs = [color_list[i]])])

    data_slice = module_copy.deepcopy(exact_data)
    data_slice.reduce([select("unit", "num_den"),
		       select("radius", samples[i][0])])
    data_slice.scale_dim("time", 1.0/60)
    g.plot(graph.data.list(data_slice.data_center_list(),
			   x = 1, y = 2,
			   title = "%s (exact)" % samples[i][1]),
	   styles = [graph.style.line(lineattrs = [color_list[i]])])

    data_slice = module_copy.deepcopy(sect_data)
    data_slice.reduce([select("unit", "num_den"),
		       select("radius", samples[i][0])])
    data_slice.scale_dim("time", 1.0/60)
    g.plot(graph.data.list(data_slice.data_center_list(),
			   x = 1, y = 2,
			   title = "%s (sect)" % samples[i][1]),
	   styles = [graph.style.symbol(symbol = graph.style.symbol.cross,
					size = 0.05,
					symbolattrs = [color_list[i]])])
    
g.writePDFfile("out/emission_history.pdf")
