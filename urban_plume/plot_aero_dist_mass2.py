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

times_hour = [1, 2, 3, 4,5, 6]
times_sec = [t * 3600 for t in times_hour]

subdir = "c"
if len(sys.argv) > 1:
    subdir = sys.argv[1]

data = pmc_var(NetCDFFile("out/urban_plume_no_coag_0001.nc"),
	       "aero", [])
data.write_summary(sys.stdout)
data.reduce([select("unit", "mass_den"),
             select("aero_species", "NO3")])

print data.data

#data.scale(1)
data.scale_dim("dry_radius", 1e6)
data.scale(1e9)

g = graph.graphxy(
    width = 10,
    x = graph.axis.log(min = 0.005,
		       max = 1,
		       title = "dry radius ($\mu$m)",
		       painter = grid_painter),
    y = graph.axis.log(min = 1.e-6,
                       max = 1.e1,
                       title = "mass density ($\mu$g/m$^3$)",
		       painter = grid_painter),
    key = graph.key.key(pos = "tr"))

for i in range(len(times_sec)):
    data_slice = module_copy.deepcopy(data)
    data_slice.reduce([select("time", times_sec[i])])
    g.plot(graph.data.list(data_slice.data_center_list(strip_zero = True),
			   x = 1, y = 2,
			   title = "%g hours" % times_hour[i]),
	   styles = [graph.style.line(lineattrs = [color_list[i], style.linewidth.THick])])

g.writePDFfile("out/t_NO3no.pdf")
