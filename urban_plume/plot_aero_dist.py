#!/usr/bin/env python
# Copyright (C) 2007 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
sys.path.append("../tool")
from pmc_data import *
from pmc_pyx import *
sys.path.append(os.path.expanduser("~/.python"))
from pyx import *

times_hour = [0, 5, 24]
times_sec = [t * 3600 for t in times_hour]

data_set = read_data_set(file_list("out/urban_plume_state", "aero.dat"),
			 [select("unit", "mass_den"),
			  sum("species")])

g = graph.graphxy(
    width = 10,
    x = graph.axis.log(title = "radius (m)",
		       painter = grid_painter),
    y = graph.axis.log(title = "mass density (kg/m$^3$)",
		       painter = grid_painter),
    key = graph.key.key(pos = "tr"))

for i in range(len(times_sec)):
    data_slice = module_copy.deepcopy(data_set)
    data_slice.reduce([select("time", times_sec[i])])
    g.plot(graph.data.list(data_slice.data_center_list(strip_zero = True),
			   x = 1, y = 2,
			   title = "%g hours" % times_hour[i]),
	   styles = [graph.style.line(lineattrs = [color_list[i]])])

g.writePDFfile("out/aero_dist.pdf")
