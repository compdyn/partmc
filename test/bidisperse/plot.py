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

# Because we have stochastic initial conditions the initial particle
# numbers may not be quite as expected. This is a problem in
# particular for the large particles, because on average we have just
# one of them, so any given run may have 0, 1, 2, etc, and anything
# other than one will not give the desired results.

small_init_num_den = module_copy.deepcopy(mc_data)
small_init_num_den.reduce([select("unit", "num_den"),
                           sum("radius", below = 5e-5),
                           select("time", 0)])
large_init_num_den = module_copy.deepcopy(mc_data)
large_init_num_den.reduce([select("unit", "num_den"),
                           sum("radius", above = 5e-5),
                           select("time", 0)])

desired_small_init_num_den = 1e9
desired_large_init_num_den = 1e5
rel_tol = 0.05
acceptable = True
if (abs((small_init_num_den.data - desired_small_init_num_den)
        / desired_small_init_num_den) > rel_tol) \
        or (abs((large_init_num_den.data - desired_large_init_num_den)
                / desired_large_init_num_den) > rel_tol):
    print ("WARNING: randomly chosen initial conditions differ by more"
           " than %g%% from the desired values, so results may be bad."
           " To fix, change the rand_init value in run_mc.spec until an"
           " acceptable run is obtained.") % (rel_tol * 100)

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

g.plot(graph.data.points(init_mc_data.data_center_list(),
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

g.plot(graph.data.points(small_mc_data.data_center_list(),
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

g.plot(graph.data.points(small_mc_data.data_center_list(strip_zero = True),
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

g.plot(graph.data.points(large_mc_data.data_center_list(),
                         x = 1, y = 2,
                         title = "Monte Carlo"),
       styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				    size = 0.05,
				    symbolattrs = [color_list[0]])])

g.writePDFfile("out/bidisperse_vol_big")

######################################################################
