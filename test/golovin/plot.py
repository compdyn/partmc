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

g_num = graph.graphxy(width = 6,
		      x = graph.axis.log(min = 1e-1,
					 max = 1e3,
					 title = "diameter ($\mu$m)",
					 painter = grid_painter),
		      y = graph.axis.log(min = 1e5,
					 max = 1e10,
					 title = "number density (\#/m$^3$)",
					 painter = grid_painter))
#		      key = graph.key.key(pos = "tr"))
g_vol = graph.graphxy(width = 6,
		      x = graph.axis.log(min = 1e-1,
					 max = 1e3,
					 title = "diameter ($\mu$m)",
					 painter = grid_painter),
		      y = graph.axis.log(min = 1e-13,
					 max = 1e-5,
					 title = "volume density (m$^3$/m$^3$)",
					 painter = grid_painter))
#		      key = graph.key.key(pos = "br"))

for i in range(len(times_sec)):
    data_slice = module_copy.deepcopy(mc_data)
    data_slice.reduce([select("unit", "num_den"),
		       select("time", times_sec[i])])
    g_num.plot(graph.data.list(data_slice.data_center_list(strip_zero = True),
			       x = 1, y = 2,
			       title = "%g mins MC" % times_min[i]),
	       styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
					    size = 0.05,
					    symbolattrs = [color.grey.black])])
    g_num.text(3.8,2.8,"0 mins",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
    g_num.text(4.2,2,"5 mins",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
    g_num.text(4.5,1.2,"10 mins",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])

    data_slice = module_copy.deepcopy(mc_data)
    data_slice.reduce([select("unit", "vol_den"),
		       select("time", times_sec[i])])
    g_vol.plot(graph.data.list(data_slice.data_center_list(strip_zero = True),
			       x = 1, y = 2,
			       title = "%g mins MC" % times_min[i]),
	       styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
					    size = 0.05,
					    symbolattrs = [color_list[i]])])
    g_vol.text(1.5,3,"0 h",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
    g_vol.text(2.5,1,"12 h",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
    g_vol.text(3.8,0.4,"24 h",[text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])

    data_slice = module_copy.deepcopy(exact_data)
    data_slice.reduce([select("unit", "num_den"),
		       select("time", times_sec[i])])
    g_num.plot(graph.data.list(data_slice.data_center_list(strip_zero = True),
			       x = 1, y = 2,
			       title = "%g mins exact" % times_min[i]),
	       styles = [graph.style.line(lineattrs = [color.grey.black])])
    
    data_slice = module_copy.deepcopy(exact_data)
    data_slice.reduce([select("unit", "vol_den"),
		       select("time", times_sec[i])])
    g_vol.plot(graph.data.list(data_slice.data_center_list(strip_zero = True),
			       x = 1, y = 2,
			       title = "%g mins exact" % times_min[i]),
	       styles = [graph.style.line(lineattrs = [color_list[i]])])
    
g_num.writePDFfile("out/golovin_num.pdf")
g_vol.writePDFfile("out/golovin_vol.pdf")

print "figure height = %.1f cm" % unit.tocm(g_num.bbox().height())
print "figure width = %.1f cm" % unit.tocm(g_num.bbox().width())
