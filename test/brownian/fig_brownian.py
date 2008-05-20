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

mc_data.scale(math.log(10.0)) # d/dln(r) to d/dlog10(r)
sect_data.scale(math.log(10.0)) # d/dln(r) to d/dlog10(r)

c = canvas.canvas()

g_vol_lin = c.insert(graph.graphxy(
    width = 6.1,
    x = graph.axis.log(min = 1e-2,
		       max = 1e0,
		       title = r"diameter ($\rm \mu m$)",
		       painter = grid_painter),
    y = graph.axis.linear(min = 0,
			  max = 1.5e-10,
			  title = r"volume density ($1$)",
			  painter = grid_painter)))

g_num_lin = c.insert(graph.graphxy(
    width = 6.1,
    ypos = g_vol_lin.height + 0.5,
    x = graph.axis.linkedaxis(g_vol_lin.axes["x"],
                              painter = graph.axis.painter.linked(gridattrs = [style.linestyle.dotted])),
    y = graph.axis.linear(min = 0,
			  max = 2e11,
			  title = r"number density ($\rm m^{-3}$)",
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

(x, y) = g_vol_lin.pos(0.075, 4e-11)
g_vol_lin.text(x, y, "0 hours", [text.halign.right, text.valign.bottom])
(x, y) = g_vol_lin.pos(0.15, 8.5e-11)
g_vol_lin.text(x, y, "12 hours", [text.halign.right, text.valign.bottom])
(x, y) = g_vol_lin.pos(0.21, 13e-11)
g_vol_lin.text(x, y, "24 hours", [text.halign.right, text.valign.bottom])

(x, y) = g_num_lin.pos(0.075, 1.6e11)
g_num_lin.text(x, y, "0 hours", [text.halign.left, text.valign.bottom])
(x, y) = g_num_lin.pos(0.18, 4e10)
g_num_lin.text(x, y, "12 hours", [text.halign.left, text.valign.bottom])
(x, y) = g_num_lin.pos(0.28, 2e10)
g_num_lin.text(x, y, "24 hours", [text.halign.left, text.valign.bottom])

c.writePDFfile("out/brownian.pdf")
print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
