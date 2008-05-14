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

time_hour = 24
composition_lower = [0,  2, 80]
composition_upper = [2, 10, 90]

withcoag_subdir = "withcoag_dry"
nocoag_subdir = "nocoag_dry"
if len(sys.argv) > 1:
    withcoag_subdir = sys.argv[1]
    nocoag_subdir = sys.argv[1]

data_no = pmc_var(NetCDFFile("out/%s/urban_plume_0001.nc" % nocoag_subdir),
	       "comp_bc",
	       [])
data_wc = pmc_var(NetCDFFile("out/%s/urban_plume_0001.nc" % withcoag_subdir),
	       "comp_bc",
	       [])

#data_no.write_summary(sys.stdout)
#data_wc.write_summary(sys.stdout)

data_no.scale_dim("composition_bc", 100)
data_wc.scale_dim("composition_bc", 100)

data_no.scale_dim("dry_radius", 2e6)
data_wc.scale_dim("dry_radius", 2e6)

data_no.scale_dim("time", 1.0/3600)
data_wc.scale_dim("time", 1.0/3600)

data_no.reduce([select("unit", "mass_den"),
                select("time", time_hour),
		sum("aero_species", without = ["H2O"])])
data_wc.reduce([select("unit", "mass_den"),
                select("time", time_hour),
		sum("aero_species", without = ["H2O"])])

#total_mass = module_copy.deepcopy(data_no)
#total_mass.reduce([sum("dry_radius"), sum("composition_bc")])
#total_mass.scale(1e9)
#print total_mass.data

data_no.scale(1e9) # kg/m^3 to ug/m^3
data_wc.scale(1e9) # kg/m^3 to ug/m^3
data_no.scale(math.log(10.0)) # d/dln(r) to d/dlog10(r)
data_wc.scale(math.log(10.0)) # d/dln(r) to d/dlog10(r)

g = graph.graphxy(
	width = 6.5,
	x = graph.axis.log(title = r'dry diameter ($\mu$m)',
                           min = 0.01, max = 2,
			   painter = grid_painter),
	y = graph.axis.log(title = r'mass density ($\rm \mu g\, m^{-3}$)',
                           min = 1e-6, max =2e2,
			   painter = grid_painter),
        key = graph.key.key(pos = "br"))

for i in range(len(composition_lower)):
    data_no_slice = module_copy.deepcopy(data_no)
    data_wc_slice = module_copy.deepcopy(data_wc)
        
    reducer1 = sum("composition_bc", above = composition_lower[i], below = composition_upper[i])
    reducers = [reducer1]
    data_no_slice.reduce(reducers)
    data_wc_slice.reduce(reducers)

    g.plot(graph.data.points(data_no_slice.data_center_list(strip_zero = True),
			   x = 1, y = 2, 
                           title = "no %d - %d" % (composition_lower[i], composition_upper[i])),
	   styles = [graph.style.line(lineattrs = [line_style_list[i],style.linewidth.Thick])])

    g.plot(graph.data.points(data_wc_slice.data_center_list(strip_zero = True),
			   x = 1, y = 2, 
                           title = "wc %d - %d" % (composition_lower[i], composition_upper[i])),
           styles = [graph.style.line(lineattrs = [line_style_list[i],style.linewidth.THick])])


g.writePDFfile("figs/bc_mixing.pdf")
print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
