#!/usr/bin/env python
# Copyright (C) 2007, 2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *

times_hour = [1, 6, 24]
times_sec = [t * 3600 for t in times_hour]

#subdir = "withcoag_dry"
if len(sys.argv) > 1:
    subdir = sys.argv[1]

data1 = pmc_var(NetCDFFile("out/withcoag_dry/urban_plume_0001.nc"),
	       "aero", [])
data1.write_summary(sys.stdout)
data1.reduce([select("unit", "mass_den"),
             select("aero_species", "BC")])

print data1.data

data1.scale_dim("dry_radius", 2e6)
data1.scale(2.303*1e9)

data2 = pmc_var(NetCDFFile("out/nocoag_dry/urban_plume_0001.nc"),
	       "aero", [])
data2.write_summary(sys.stdout)
data2.reduce([select("unit", "mass_den"),
             select("aero_species", "BC")])

print data2.data

data2.scale_dim("dry_radius", 2e6)
data2.scale(2.303*1e9)

g = graph.graphxy(
    width = 10,
    x = graph.axis.log(min = 0.01,
		       max = 2,
		       title = "dry diameter ($\mu$m)",
		       painter = grid_painter),
    y = graph.axis.log(min = 1.e-5,
                       max = 1.e1,
                       title = "mass density ($\mu$g/m$^3$)",
		       painter = grid_painter),
    key = graph.key.key(pos = "tl"))

for i in range(len(times_sec)):
    data1_slice = module_copy.deepcopy(data1)
    data1_slice.reduce([select("time", times_sec[i])])
    g.plot(graph.data.list(data1_slice.data_center_list(strip_zero = True),
			   x = 1, y = 2,
			   title = "%g hours" % times_hour[i]),
	   styles = [graph.style.line(lineattrs = [color_list[i], style.linewidth.THick])])
data2_slice = module_copy.deepcopy(data2)
data2_slice.reduce([select("time", times_sec[2])])
g.plot(graph.data.list(data2_slice.data_center_list(strip_zero = True),
			   x = 1, y = 2,
			   title = "%g hours, no coag" % times_hour[2]),
	   styles = [graph.style.line(lineattrs = [color_list[4], style.linewidth.THick])])

g.writePDFfile("out/withcoag_dry/aero_dist_BC.pdf")
