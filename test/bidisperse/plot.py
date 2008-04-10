#!/usr/bin/env python
# Copyright (C) 2007, 2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../../tool")
from pmc_data_nc import *
from pmc_pyx import *

mc_data = pmc_var(NetCDFFile("out/bidisperse_mc_0001.nc"),
		  "aero",
		  [sum("aero_species")])

init_mc_data = module_copy.deepcopy(mc_data)
init_mc_data.reduce([select("unit", "vol_den"),
		     select("time", 0)])

small_mc_data = module_copy.deepcopy(mc_data)
small_mc_data.reduce([select("unit", "num_den"),
		      sum("radius", below = 5e-5)])

large_mc_data = module_copy.deepcopy(mc_data)
large_mc_data.reduce([select("unit", "vol_den"),
		      sum("radius", above = 5e-5)])

######################################################################

g = graph.graphxy(
    width = 10,
    x = graph.axis.log(
	min = 1e-6,
	max = 1e-3,
	title = "radius (m)",
	painter = grid_painter),
    y = graph.axis.linear(
	title = "volume density (m$^3$/m$^3$)",
	painter = grid_painter),
    key = graph.key.key(pos = "tr"))

g.plot(graph.data.list(init_mc_data.data_center_list(),
		       x = 1, y = 2,
		       title = "initial distribution"),
       styles = [graph.style.line(lineattrs = [color_list[0]])])

g.writePDFfile("out/bidisperse_init_dist")

######################################################################

g = graph.graphxy(
    width = 10,
    x = graph.axis.linear(
	title = "time (s)",
	painter = grid_painter),
    y = graph.axis.linear(
	title = "number density of small particles (\#/m$^3$)",
	painter = grid_painter),
    key = graph.key.key(pos = "tr"))

g.plot(graph.data.file("out/bidisperse_ode_counts.d",
		       x = 1, y = 2,
		       title = "ODE"),
       styles = [graph.style.line(lineattrs = [color_list[1]])])

g.plot(graph.data.list(small_mc_data.data_center_list(),
		       x = 1, y = 2,
		       title = "Monte Carlo"),
       styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				    size = 0.05,
				    symbolattrs = [color_list[0]])])

g.writePDFfile("out/bidisperse_num_small")

######################################################################

g = graph.graphxy(
    width = 10,
    x = graph.axis.linear(
	title = "time (s)",
	painter = grid_painter),
    y = graph.axis.log(
	title = "number density of small particles (\#/m$^3$)",
	painter = grid_painter),
    key = graph.key.key(pos = "tr"))

g.plot(graph.data.file("out/bidisperse_ode_counts.d",
		       x = 1, y = 2,
		       title = "ODE"),
       styles = [graph.style.line(lineattrs = [color_list[1]])])

g.plot(graph.data.list(small_mc_data.data_center_list(strip_zero = True),
		       x = 1, y = 2,
		       title = "Monte Carlo"),
       styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				    size = 0.05,
				    symbolattrs = [color_list[0]])])

g.writePDFfile("out/bidisperse_num_small_log")

######################################################################

g = graph.graphxy(
    width = 10,
    x = graph.axis.linear(
	title = "time (s)",
	painter = grid_painter),
    y = graph.axis.linear(
	title = "volume density of big particles (m$^3$/m$^3$)",
	painter = grid_painter),
    key = graph.key.key(pos = "br"))

g.plot(graph.data.file("out/bidisperse_ode_counts.d",
		       x = 1, y = 3,
		       title = "ODE"),
       styles = [graph.style.line(lineattrs = [color_list[1]])])

g.plot(graph.data.list(large_mc_data.data_center_list(),
		       x = 1, y = 2,
		       title = "Monte Carlo"),
       styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				    size = 0.05,
				    symbolattrs = [color_list[0]])])

g.writePDFfile("out/bidisperse_vol_big")

######################################################################
