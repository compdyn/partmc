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

mc_data = pmc_var(NetCDFFile("out/brown_mc_0001.nc"),
		  "aero",
		  [sum("aero_species")])
sect_data = pmc_var(NetCDFFile("out/brown_sect_0001.nc"),
		    "aero",
		    [sum("aero_species")])

mc_data.scale_dim("radius", 2e6)
sect_data.scale_dim("radius",2e6)

c = canvas.canvas()

g_vol_lin = c.insert(graph.graphxy(
    width = 6,
    x = graph.axis.log(min = 1e-2,
		       max = 1e0,
		       title = "diameter ($\mu$m)",
		       painter = grid_painter),
    y = graph.axis.linear(min = 0,
			  max = 4e-11,
			  title = "volume density (m$^3$/m$^3$)",
			  painter = grid_painter)))

g_num_lin = c.insert(graph.graphxy(
    width = 6,
    ypos = g_vol_lin.height + 0.5,
    x = graph.axis.linkedaxis(g_vol_lin.axes["x"],
                              painter = graph.axis.painter.linked(gridattrs = [style.linestyle.dotted])),
    y = graph.axis.linear(min = 0,
			  max = 5e10,
			  title = "number density (\#/m$^3$)",
			  painter = grid_painter)))

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

    data_slice = module_copy.deepcopy(sect_data)
    data_slice.reduce([select("unit", "num_den"),
		       select("time", times_sec[i])])

    g_num_lin.plot(
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

g_num_lin.text(g_num_lin.xpos + 1.5,
               g_num_lin.ypos + 2.9,
               "0 h",
               [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g_num_lin.text(g_num_lin.xpos + 2.5,
               g_num_lin.ypos + 1.1,
               "12 h",
               [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g_num_lin.text(g_num_lin.xpos + 3,
               g_num_lin.ypos + 0.1,
               "24 h",
               [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])

g_vol_lin.text(g_vol_lin.xpos + 2.3,
               g_vol_lin.ypos + 1.1,
               "0 h",
               [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g_vol_lin.text(g_vol_lin.xpos + 2.8,
               g_vol_lin.ypos + 2.2,
               "12 h",
               [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g_vol_lin.text(g_vol_lin.xpos + 3.3,
               g_vol_lin.ypos + 3.3,
               "24 h",
               [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
#g_vol_lin.stroke(path.line(2.2, 0.6, 2.95, 0.5))

c.writePDFfile("out/brownian.pdf")
print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
