#!/usr/bin/env python
# Copyright (C) 2007 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
sys.path.append("../../tool")
from pmc_data_nc import *
from pmc_pyx import *
sys.path.append(os.path.expanduser("~/.python"))
from pyx import *
from Scientific.IO.NetCDF import *

times_min = [0, 5, 10]
times_sec = [t * 60 for t in times_min]

mc_data = pmc_var(NetCDFFile("out/sedi_mc_state_0001.nc"),
		  "aero",
		  [sum("aero_species")])
sect_data = pmc_var(NetCDFFile("out/sedi_sect_0001.nc"),
		    "aero",
		    [sum("aero_species")])

g_num_lin = graph.graphxy(
    width = 10,
    x = graph.axis.log(min = 1e-7,
		       max = 1e-2,
		       title = "radius (m)",
		       painter = grid_painter),
    y = graph.axis.linear(min = 0,
			  max = 1.25e9,
			  title = "number density (\#/m$^3$)",
			  painter = grid_painter),
    key = graph.key.key(pos = "tr"))
g_num_log = graph.graphxy(
    width = 10,
    x = graph.axis.log(min = 1e-7,
		       max = 1e-2,
		       title = "radius (m)",
		       painter = grid_painter),
    y = graph.axis.log(min = 1e2,
		       max = 1e10,
		       title = "number density (\#/m$^3$)",
		       painter = grid_painter),
    key = graph.key.key(pos = "tr"))
g_vol_lin = graph.graphxy(
    width = 10,
    x = graph.axis.log(min = 1e-7,
		       max = 1e-2,
		       title = "radius (m)",
		       painter = grid_painter),
    y = graph.axis.linear(min = 0,
			  max = 8e-6,
			  title = "volume density (m$^3$/m$^3$)",
			  painter = grid_painter),
    key = graph.key.key(pos = "tl"))
g_vol_log = graph.graphxy(
    width = 10,
    x = graph.axis.log(min = 1e-7,
		       max = 1e-2,
		       title = "radius (m)",
		       painter = grid_painter),
    y = graph.axis.log(min = 1e-20,
		       max = 1e-3,
		       title = "volume density (m$^3$/m$^3$)",
		       painter = grid_painter),
    key = graph.key.key(pos = "bl", hdist = 1.8))

for i in range(len(times_sec)):
    data_slice = module_copy.deepcopy(mc_data)
    data_slice.reduce([select("unit", "num_den"),
		       select("time", times_sec[i])])
    g_num_lin.plot(
	graph.data.list(data_slice.data_center_list(strip_zero = True),
			x = 1, y = 2,
			title = "%g mins MC" % times_min[i]),
	styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				     size = 0.05,
				     symbolattrs = [color_list[i]])])
    g_num_log.plot(
	graph.data.list(data_slice.data_center_list(strip_zero = True),
			x = 1, y = 2,
			title = "%g mins MC" % times_min[i]),
	styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				     size = 0.05,
				     symbolattrs = [color_list[i]])])

    data_slice = module_copy.deepcopy(mc_data)
    data_slice.reduce([select("unit", "vol_den"),
		       select("time", times_sec[i])])
    g_vol_lin.plot(
	graph.data.list(data_slice.data_center_list(strip_zero = True),
			x = 1, y = 2,
			title = "%g mins MC" % times_min[i]),
	styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				     size = 0.05,
				     symbolattrs = [color_list[i]])])
    g_vol_log.plot(
	graph.data.list(data_slice.data_center_list(strip_zero = True),
			x = 1, y = 2,
			title = "%g mins MC" % times_min[i]),
	styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				     size = 0.05,
				     symbolattrs = [color_list[i]])])

    data_slice = module_copy.deepcopy(sect_data)
    data_slice.reduce([select("unit", "num_den"),
		       select("time", times_sec[i])])
    g_num_lin.plot(
	graph.data.list(data_slice.data_center_list(strip_zero = True),
			x = 1, y = 2,
			title = "%g mins sect" % times_min[i]),
	styles = [graph.style.line(lineattrs = [color_list[i]])])
    g_num_log.plot(
	graph.data.list(data_slice.data_center_list(strip_zero = True),
			x = 1, y = 2,
			title = "%g mins sect" % times_min[i]),
	styles = [graph.style.line(lineattrs = [color_list[i]])])
    
    data_slice = module_copy.deepcopy(sect_data)
    data_slice.reduce([select("unit", "vol_den"),
		       select("time", times_sec[i])])
    g_vol_lin.plot(
	graph.data.list(data_slice.data_center_list(strip_zero = True),
			x = 1, y = 2,
			title = "%g mins sect" % times_min[i]),
	styles = [graph.style.line(lineattrs = [color_list[i]])])
    g_vol_log.plot(
	graph.data.list(data_slice.data_center_list(strip_zero = True),
			x = 1, y = 2,
			title = "%g mins sect" % times_min[i]),
	styles = [graph.style.line(lineattrs = [color_list[i]])])
    
g_num_lin.writePDFfile("out/sedi_num_lin.pdf")
g_num_log.writePDFfile("out/sedi_num_log.pdf")
g_vol_lin.writePDFfile("out/sedi_vol_lin.pdf")
g_vol_log.writePDFfile("out/sedi_vol_log.pdf")
