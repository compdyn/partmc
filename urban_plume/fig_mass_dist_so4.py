#!/usr/bin/env python
# Copyright (C) 2007, 2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *

times_hour = [1, 6, 24]
times_sec = [t * 3600 for t in times_hour]
line_style_order = [2, 1, 0]

data1 = pmc_var(NetCDFFile("out/urban_plume_with_coag_0001.nc"),
	       "aero", [])
#data1.write_summary(sys.stdout)
data1.reduce([select("unit", "mass_den"),
             select("aero_species", "SO4")])

data1.scale_dim("dry_radius", 1e6) # m to um
data1.scale_dim("dry_radius", 2.0) # radius to diameter
data1.scale(1e9) # kg/m^3 to ug/m^3
data1.scale(math.log(10.0)) # d/dln(r) to d/dlog10(r)

data2 = pmc_var(NetCDFFile("out/urban_plume_no_coag_0001.nc"),
	       "aero", [])
#data2.write_summary(sys.stdout)
data2.reduce([select("unit", "mass_den"),
             select("aero_species", "SO4")])

data2.scale_dim("dry_radius", 1e6) # m to um
data2.scale_dim("dry_radius", 2.0) # radius to diameter
data2.scale(1e9) # kg/m^3 to ug/m^3
data2.scale(math.log(10.0)) # d/dln(r) to d/dlog10(r)

g = graph.graphxy(
    width = 6.4,
    x = graph.axis.log(min = 0.01,
		       max = 2,
		       title = "dry diameter ($\mu$m)",
		       painter = major_grid_painter),
    y = graph.axis.log(min = 1.e-6,
                       max = 2.e1,
                       title = "mass density ($\mu$g/m$^3$)",
		       painter = major_grid_painter),
    key = graph.key.key(pos = "bc"))

for i in range(len(times_sec)):
    data1_slice = module_copy.deepcopy(data1)
    data1_slice.reduce([select("time", times_sec[i])])
    g.plot(graph.data.list(data1_slice.data_center_list(strip_zero = True),
			   x = 1, y = 2,
			   title = "%g hours" % times_hour[i]),
	   styles = [graph.style.line(lineattrs = [line_style_list[line_style_order[i]], style.linewidth.Thick])])
data2_slice = module_copy.deepcopy(data2)
data2_slice.reduce([select("time", times_sec[2])])
g.plot(graph.data.list(data2_slice.data_center_list(strip_zero = True),
			   x = 1, y = 2,
			   title = "%g hours, no coag" % times_hour[2]),
	   styles = [graph.style.line(lineattrs = [line_style_list[0], style.linewidth.THick])])

g.writePDFfile("figs/mass_dist_so4.pdf")
print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
