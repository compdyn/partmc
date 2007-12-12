#!/usr/bin/env python
# Copyright (C) 2007 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../../tool")
from pmc_data_nc import *
from pmc_pyx import *

times_hour = [1, 6, 12, 18, 24]

data = pmc_var(NetCDFFile("out/urban_plume_state_0001.nc"),
	       "optic_scatter",
	       [])
data.write_summary(sys.stdout)

data.reduce([select("unit", "num_den"),
		 sum("aero_species")])
data.scale_dim("radius", 1e6)
data.scale_dim("time", 1.0/3600)
data.scale_dim("scatter_cross_section_area", 1e12)

for i in range(len(times_hour)):
    g = graph.graphxy(
	width = 10,
	x = graph.axis.log(title = r'radius ($\mu$m)',
			   painter = grid_painter),
	y = graph.axis.log(title = r'scattering cross sectional area ($\mu {\rm m}^2$)',
			   painter = grid_painter))
    data_slice = module_copy.deepcopy(data)
    data_slice.reduce([select("time", times_hour[i])])
    g.plot(graph.data.list(data_slice.data_2d_list(strip_zero = True),
			   xmin = 1, xmax = 2, ymin = 3, ymax = 4, color = 5),
	   styles = [graph.style.rect(rainbow_palette)])
    g.writePDFfile("out/aero_optic_scatter_%d.pdf" % times_hour[i])
