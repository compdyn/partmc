#!/usr/bin/env python
# Copyright (C) 2007-2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *

time_hour = 24

data = pmc_var(NetCDFFile("out/urban_plume_with_coag_0001.nc"),
	       "n_orig",
	       [])

data.reduce([select("unit", "num_den"),
		 sum("aero_species")])
data.scale_dim("dry_radius", 1e6)
data.scale_dim("time", 1.0/3600)
data.scale(math.log(10)) # d/dln(r) to d/dlog10(r)

data.dim_by_name("n_orig_part").grid_centers = data.dim_by_name("n_orig_part").grid_centers - 1
data.dim_by_name("n_orig_part").grid_edges = data.dim_by_name("n_orig_part").grid_edges - 1

g = graph.graphxy(
    width = 6.9,
    x = graph.axis.log(min = 1e-3,
                       max = 1e+0,
                       title = r'dry radius ($\rm \mu m$)',
                       painter = grid_painter),
    y = graph.axis.linear(title = 'number of coagulation events',
                          painter = grid_painter))
data_slice = module_copy.deepcopy(data)
data_slice.reduce([select("time", time_hour)])
#    min_val = data_slice.data.min()
#    max_val = data_slice.data.max()
min_val = 0
max_val = 1e9
g.plot(graph.data.points(data_slice.data_2d_list(strip_zero = True,
                                               min = min_val,
                                               max = max_val),
                       xmin = 1, xmax = 2, ymin = 3,
                       ymax = 4, color = 5),
       styles = [graph.style.rect(gray_palette)])
add_horiz_color_bar(g,
              min = min_val,
              max = max_val,
              title = r"number density ($\rm m^{-3}$)",
              palette = gray_palette,
              bar_offset = 0.6)
g.writePDFfile("figs/n_orig.pdf")
print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
