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

times_min = [0, 720, 1440]
times_sec = [t * 60 for t in times_min]
times_hour = [t / 60  for t in times_min]

mc_data = pmc_var(NetCDFFile("out/brown_mc_avg.nc"),
		  "aero",
		  [sum("aero_species")])
sect_data = pmc_var(NetCDFFile("out/brown_sect_0001.nc"),
		    "aero",
		    [sum("aero_species")])

mc_data.scale_dim("radius", 2e6)
sect_data.scale_dim("radius",2e6)

g_num_lin = graph.graphxy(
    width = 6,
    x = graph.axis.log(min = 1e-2,
		       max = 1e0,
		       title = "diameter ($\mu$m)",
		       painter = grid_painter),
    y = graph.axis.linear(min = 0,
			  max = 8e10,
			  title = "number density (\#/m$^3$)",
			  painter = grid_painter))
g_num_log = graph.graphxy(
    width = 6,
    x = graph.axis.log(min = 1e-2,
		       max = 1e0,
		       title = "diameter ($\mu$m)",
		       painter = grid_painter),
    y = graph.axis.log(min = 1e6,
		       max = 1e12,
		       title = "number density (\#/m$^3$)",
		       painter = major_grid_painter))
g_vol_lin = graph.graphxy(
    width = 6,
    x = graph.axis.log(min = 1e-2,
		       max = 1e0,
		       title = "diameter ($\mu$m)",
		       painter = grid_painter),
    y = graph.axis.linear(min = 0,
			  max = 6e-11,
			  title = "volume density (m$^3$/m$^3$)",
			  painter = grid_painter))
g_vol_log = graph.graphxy(
    width = 6,
    x = graph.axis.log(min = 1e-2,
		       max = 1e0,
		       title = "diameter ($\mu$m)",
		       painter = grid_painter),
    y = graph.axis.log(min = 1e-17,
		       max = 1e-10,
		       title = "volume density (m$^3$/m$^3$)",
		       painter = grid_painter))

for i in range(len(times_sec)):
    data_slice = module_copy.deepcopy(mc_data)
    data_slice.reduce([select("unit", "num_den"),
		       select("time", times_sec[i])])
    g_num_lin.plot(
	graph.data.points(data_slice.data_center_list(strip_zero = True),
                          x = 1, y = 2,
                          title = "%g hours MC" % times_hour[i]),
	styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				     size = 0.05,
				     symbolattrs = [color.grey.black])])

    g_num_log.plot(
	graph.data.points(data_slice.data_center_list(strip_zero = True),
                          x = 1, y = 2,
                          title = "%g hours MC" % times_hour[i]),
	styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				     size = 0.05,
				     symbolattrs = [color.grey.black])])

    data_slice = module_copy.deepcopy(mc_data)
    data_slice.reduce([select("unit", "vol_den"),
		       select("time", times_sec[i])])
    g_vol_lin.plot(
	graph.data.points(data_slice.data_center_list(strip_zero = True),
                          x = 1, y = 2,
                          title = "%g hours MC" % times_hour[i]),
	styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				     size = 0.05,
				     symbolattrs = [color.grey.black])])

    g_vol_log.plot(
	graph.data.points(data_slice.data_center_list(strip_zero = True),
                          x = 1, y = 2,
                          title = "%g hours MC" % times_hour[i]),
	styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				     size = 0.05,
				     symbolattrs = [color.grey.black])])

    data_slice = module_copy.deepcopy(sect_data)
    data_slice.reduce([select("unit", "num_den"),
		       select("time", times_sec[i])])

    g_num_lin.plot(
	graph.data.points(data_slice.data_center_list(strip_zero = True),
                          x = 1, y = 2,
                          title = "%g hours sect" % times_hour[i]),
	styles = [graph.style.line(lineattrs = [color.grey.black])])
    g_num_log.plot(
	graph.data.points(data_slice.data_center_list(strip_zero = True),
                          x = 1, y = 2,
                          title = "%g hours sect" % times_hour[i]),
	styles = [graph.style.line(lineattrs = [color.grey.black])])
    
    data_slice = module_copy.deepcopy(sect_data)
    data_slice.reduce([select("unit", "vol_den"),
		       select("time", times_sec[i])])
    g_vol_lin.plot(
	graph.data.points(data_slice.data_center_list(strip_zero = True),
                          x = 1, y = 2,
                          title = "%g hours sect" % times_hour[i]),
	styles = [graph.style.line(lineattrs = [color.grey.black])])
    g_vol_log.plot(
	graph.data.points(data_slice.data_center_list(strip_zero = True),
                          x = 1, y = 2,
                          title = "%g hours sect" % times_hour[i]),
	styles = [graph.style.line(lineattrs = [color.grey.black])])
    
g_num_lin.text(2.7,3,"0 h",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g_num_lin.text(3.8,0.8,"12 h",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g_num_lin.text(4.3,0.4,"24 h",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g_num_log.text(0.5,2.8,"0 h",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g_num_log.text(1.2,2,"12 h",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g_num_log.text(1.9,0.6,"24 h",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g_vol_lin.text(2.1,1.1,"0 h",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g_vol_lin.text(2.9,2.4,"12 h",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g_vol_lin.text(3.3,3.3,"24 h",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g_vol_log.text(1.1,2.8,"0 h",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g_vol_log.text(0.8,1.3,"12 h",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g_vol_log.text(1.8,0.4,"24 h",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])

g_num_lin.writePDFfile("out/brown_num_lin.pdf")
g_num_log.writePDFfile("out/brown_num_log.pdf")
g_vol_lin.writePDFfile("out/brown_vol_lin.pdf")
g_vol_log.writePDFfile("out/brown_vol_log.pdf")
