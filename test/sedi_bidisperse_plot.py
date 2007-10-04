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

mc_data = pmc_var(NetCDFFile("out/sedi_bidisperse_mc_state_0001.nc"),
		  "aero",
		  [sum("aero_species")])

init_mc_data = module_copy.deepcopy(mc_data)
init_mc_data.reduce([select("unit", "vol_den"),
		     select("time", 0)])

small_mc_data = module_copy.deepcopy(mc_data)
small_mc_data.reduce([select("unit", "num_den"),
		      sum_below("radius", 5e-5)])

large_mc_data = module_copy.deepcopy(mc_data)
large_mc_data.reduce([select("unit", "vol_den"),
		      sum_above("radius", 5e-5)])

######################################################################

g = graph.graphxy(
    width = 10,
    x = graph.axis.log(title = "radius (m)",
		       painter = grid_painter),
    y = graph.axis.linear(title = "volume density (m$^3$/m$^3$)",
			  painter = grid_painter),
    key = graph.key.key(pos = "tr"))

g.plot(graph.data.list(init_mc_data.data_center_list(),
		       x = 1, y = 2,
		       title = "initial distribution"),
       styles = [graph.style.line(lineattrs = [color_list[0]])])
g.writePDFfile("out/sedi_bidisperse_init_dist")

######################################################################

g = graph.graphxy(width = 10,
		  x = graph.axis.linear(title = "time (s)",
					painter = grid_painter),
		  y = graph.axis.linear(title = "number density of small particles (\#/m$^3$)",
					painter = grid_painter),
		  key = graph.key.key(pos = "tr"))
g.plot(graph.data.file("out/sedi_bidisperse_mc_counts.d",
		       x = 1, y = 2,
		       title = "Monte Carlo"),
       styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				    size = 0.05,
				    symbolattrs = [color_list[0]])])
g.plot(graph.data.file("out/sedi_bidisperse_ode_counts.d",
		       x = 1, y = 2,
		       title = "ODE"),
       styles = [graph.style.line(lineattrs = [color_list[1]])])
g.plot(graph.data.list(small_mc_data.data_center_list(),
		       x = 1, y = 2,
		       title = "Monte Carlo NC"),
       styles = [graph.style.symbol(symbol = graph.style.symbol.cross,
				    size = 0.05,
				    symbolattrs = [color_list[1]])])
g.writePDFfile("out/sedi_bidisperse_num_small")

######################################################################

file = open("out/sedi_bidisperse_mc_counts.d")
data = []
for line in file:
    data.append([float(x) for x in line.split()])
data_small = [[d[0], d[1]] for d in data if d[1] > 0]
file.close()

g = graph.graphxy(width = 10,
		  x = graph.axis.linear(title = "time (s)",
					painter = grid_painter),
		  y = graph.axis.log(title = "number density of small particles (\#/m$^3$)",
				     painter = grid_painter),
		  key = graph.key.key(pos = "tr"))
g.plot(graph.data.list(data_small, x = 1, y = 2,
		       title = "Monte Carlo"),
       styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				    size = 0.05,
				    symbolattrs = [color_list[0]])])
g.plot(graph.data.file("out/sedi_bidisperse_ode_counts.d",
		       x = 1, y = 2,
		       title = "ODE"),
       styles = [graph.style.line(lineattrs = [color_list[1]])])
g.plot(graph.data.list(small_mc_data.data_center_list(strip_zero = True),
		       x = 1, y = 2,
		       title = "Monte Carlo NC"),
       styles = [graph.style.symbol(symbol = graph.style.symbol.cross,
				    size = 0.05,
				    symbolattrs = [color_list[1]])])
g.writePDFfile("out/sedi_bidisperse_num_small_log")

######################################################################

g = graph.graphxy(width = 10,
		  x = graph.axis.linear(title = "time (s)",
					painter = grid_painter),
		  y = graph.axis.linear(title = "volume density of big particles (m$^3$/m$^3$)",
					painter = grid_painter),
		  key = graph.key.key(pos = "br"))
g.plot(graph.data.file("out/sedi_bidisperse_mc_counts.d",
		       x = 1, y = 3,
		       title = "Monte Carlo"),
       styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				    size = 0.05,
				    symbolattrs = [color_list[0]])])
g.plot(graph.data.file("out/sedi_bidisperse_ode_counts.d",
		       x = 1, y = 3,
		       title = "ODE"),
       styles = [graph.style.line(lineattrs = [color_list[1]])])
g.plot(graph.data.list(large_mc_data.data_center_list(),
		       x = 1, y = 2,
		       title = "Monte Carlo NC"),
       styles = [graph.style.symbol(symbol = graph.style.symbol.cross,
				    size = 0.05,
				    symbolattrs = [color_list[1]])])
g.writePDFfile("out/sedi_bidisperse_vol_big")
