#!/usr/bin/env python
# Copyright (C) 2007 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *

times_hour = [6, 12, 18, 24]
times_sec = [t * 3600 for t in times_hour]
#times_sec = [ 0, 300, 600, 900]
data1 = pmc_var(NetCDFFile("out/testcase_withcoag/urban_plume_state_0001.nc"),
	       "aero",
	       [select("unit", "num_den"),
		sum("aero_species")])
data1.write_summary(sys.stdout)

data1.scale(1)
data1.scale_dim("radius", 1e6)

data2 = pmc_var(NetCDFFile("out/testcase_nocoag/urban_plume_state_0001.nc"),
               "aero",
               [select("unit", "num_den"),
                sum("aero_species")])
data2.write_summary(sys.stdout)

data2.scale(1)
data2.scale_dim("radius", 1e6)


g = graph.graphxy(
    width = 10,
    x = graph.axis.log(min = 0.01,
		       max = 2,
		       title = "radius ($\mu$m)",
		       painter = grid_painter),
    y = graph.axis.log(min = 1e7,
                       max = 2e10,
                       title = "number density (1/m$^3$)",
		       painter = grid_painter),
    key = graph.key.key(pos = "tr"))

for i in range(len(times_sec)):
    data1_slice = module_copy.deepcopy(data1)
    data1_slice.reduce([select("time", times_sec[i])])
    g.plot(graph.data.list(data1_slice.data_center_list(strip_zero = True),
			   x = 1, y = 2,
			   title = "%g hours with" % times_hour[i]),
	   styles = [graph.style.line(lineattrs = [color_list[i], style.linewidth.THick])])
    data2_slice = module_copy.deepcopy(data2)
    data2_slice.reduce([select("time", times_sec[i])])
    g.plot(graph.data.list(data2_slice.data_center_list(strip_zero = True),
                           x = 1, y = 2,
                           title = "%g hours no" % times_hour[i]),
           styles = [graph.style.line(lineattrs = [color_list[i]])])

g.writePDFfile("out/testcase_withcoag/aero_dist_num_ww.pdf")
