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

env = ["H"]

subdir = "."
if len(sys.argv) > 1:
    subdir = sys.argv[1]

data = pmc_var(NetCDFFile("out/%s/urban_plume_0001.nc" % subdir),
	       "env_state",
	       [])

data.write_summary(sys.stdout)

data.reduce([select("env", "height")])
print data.data
g = graph.graphxy(
    width = 10,
    x = graph.axis.linear(title = "time (hour)",
			  painter = grid_painter),
    y = graph.axis.linear(title = "mixing height (m)",
			  painter = grid_painter),
    key = graph.key.key(pos = "tr"))

data_slice = module_copy.deepcopy(data)
print data_slice.data
data_slice.scale_dim("time", 1.0/3600)
#data_slice.scale(100.)
g.plot(graph.data.list(data_slice.data_center_list(),
			   x = 1, y = 2,
                           title = "H"),
             styles = [graph.style.line(lineattrs = [color_list[1]])])

g.writePDFfile("out/%s/env_h.pdf" % subdir)
