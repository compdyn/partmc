#!/usr/bin/env python
# Copyright (C) 2007-2008 Matthew West, Nicole Riemer
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, glob
import copy as module_copy
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
sys.path.append(os.path.expanduser("~/.python"))
from pyx import *
from Scientific.IO.NetCDF import *

times_hour = [1, 6, 12, 18, 24]

subdir = "."
if len(sys.argv) > 1:
    subdir = sys.argv[1]

data = pmc_var(NetCDFFile("out/%s/urban_plume_0001.nc" % subdir),
	       "comp_bc",
	       [])
data.write_summary(sys.stdout)
data.scale_dim("composition", 100)
data.scale_dim("radius", 1e6)
data.scale_dim("time", 1.0/3600)

data.reduce([select("unit", "num_den"),
		 sum("aero_species"),
                 sum("composition", above = 0, below = 20)])

g = graph.graphxy(
	width = 10,
	x = graph.axis.log(title = r'radius ($\mu$m)',
                           min = 0.01, max = 1,
			   painter = grid_painter),
	y = graph.axis.log(min = 1.e6, max = 1.e10,
			      title = r'number density ($\rm m^{-3}$)',
			      painter = grid_painter),
        key = graph.key.key(pos = "tr"))

for i in range(len(times_hour)):
    data_slice = module_copy.deepcopy(data)
    data_slice.reduce([select("time", times_hour[i])])
    g.plot(graph.data.list(data_slice.data_center_list(strip_zero = True),
			   x = 1, y = 2, 
                           title = "%g hours" % times_hour[i]),
	   styles = [graph.style.line(lineattrs = [color_list[i]])])

g.writePDFfile("out/%s/aero_comp_bc1d0-20.pdf" % subdir)
