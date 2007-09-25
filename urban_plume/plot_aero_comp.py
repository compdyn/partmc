#!/usr/bin/env python
# Copyright (C) 2007 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, glob
import copy as module_copy
sys.path.append("../tool")
from pmc_data import *
from pmc_pyx import *
sys.path.append(os.path.expanduser("~/.python"))
from pyx import *

times_hour = [1, 6, 24]

data_set = read_data_set(glob.glob("out/urban_plume_state*comp_bc.dat"),
			 [])
data_set.write_summary(sys.stdout)

data_set.reduce([select("unit", "num_den"),
		 sum("species")])
data_set.scale_dim("composition", 100)
data_set.scale_dim("radius", 1e6)
data_set.scale_dim("time", 1.0/3600)

for i in range(len(times_hour)):
    g = graph.graphxy(
	width = 10,
	x = graph.axis.log(title = r'radius ($\mu$m)',
			   painter = grid_painter),
	y = graph.axis.linear(min = 0,
			      max = 100,
			      title = 'soot volume fraction',
			      texter = graph.axis.texter.decimal(suffix
								 = r"\%"),
			      painter = grid_painter))
    data_slice = module_copy.deepcopy(data_set)
    data_slice.reduce([select("time", times_hour[i])])
    g.plot(graph.data.list(data_slice.data_2d_list(strip_zero = True,
						   flip_axes = True),
			   xmin = 1, xmax = 2, ymin = 3, ymax = 4, color = 5),
	   styles = [graph.style.rect(rainbow_palette)])
    g.writePDFfile("out/aero_comp_bc_%d.pdf" % times_hour[i])
