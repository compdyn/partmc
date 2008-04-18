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

times_hour = [6]

data1 = pmc_var(NetCDFFile("out/testcase_withcoag/urban_plume_state_0001.nc"),
	       "comp_bc",
	       [])
data1.write_summary(sys.stdout)
data1.scale_dim("composition", 100)
data1.scale_dim("radius", 1e6)
data1.scale_dim("time", 1.0/3600)

data1.reduce([select("unit", "num_den"),
		 sum("aero_species"),
                 sum("composition", above = 0, below = 20)])

data2 = pmc_var(NetCDFFile("out/testcase_withcoag/urban_plume_state_0001.nc"),
               "comp_bc",
               [])
data2.write_summary(sys.stdout)
data2.scale_dim("composition", 100)
data2.scale_dim("radius", 1e6)
data2.scale_dim("time", 1.0/3600)

data2.reduce([select("unit", "num_den"),
                 sum("aero_species"),
                 sum("composition", above = 20, below = 40)])

data3 = pmc_var(NetCDFFile("out/testcase_withcoag/urban_plume_state_0001.nc"),
               "comp_bc",
               [])
data3.write_summary(sys.stdout)
data3.scale_dim("composition", 100)
data3.scale_dim("radius", 1e6)
data3.scale_dim("time", 1.0/3600)

data3.reduce([select("unit", "num_den"),
                 sum("aero_species"),
                 sum("composition", above = 40, below = 60)])

data4 = pmc_var(NetCDFFile("out/testcase_withcoag/urban_plume_state_0001.nc"),
               "comp_bc",
               [])
data4.write_summary(sys.stdout)
data4.scale_dim("composition", 100)
data4.scale_dim("radius", 1e6)
data4.scale_dim("time", 1.0/3600)

data4.reduce([select("unit", "num_den"),
                 sum("aero_species"),
                 sum("composition", above = 60, below = 80)])

data5 = pmc_var(NetCDFFile("out/testcase_withcoag/urban_plume_state_0001.nc"),
               "comp_bc",
               [])
data5.write_summary(sys.stdout)
data5.scale_dim("composition", 100)
data5.scale_dim("radius", 1e6)
data5.scale_dim("time", 1.0/3600)

data5.reduce([select("unit", "num_den"),
                 sum("aero_species"),
                 sum("composition", above = 80, below = 100)])


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
    data1_slice = module_copy.deepcopy(data1)
    data1_slice.reduce([select("time", times_hour[i])])
    g.plot(graph.data.list(data1_slice.data_center_list(strip_zero = True),
			   x = 1, y = 2, 
                           title = "0-20"),
	   styles = [graph.style.line(lineattrs = [color_list[i+1]])])
    data2_slice = module_copy.deepcopy(data2)
    data2_slice.reduce([select("time", times_hour[i])])
    g.plot(graph.data.list(data2_slice.data_center_list(strip_zero = True),
                           x = 1, y = 2,
                           title = "20-40"),
           styles = [graph.style.line(lineattrs = [color_list[i+2]])])
    data3_slice = module_copy.deepcopy(data3)
    data3_slice.reduce([select("time", times_hour[i])])
    g.plot(graph.data.list(data3_slice.data_center_list(strip_zero = True),
                           x = 1, y = 2,
                           title = "40-60"),
           styles = [graph.style.line(lineattrs = [color_list[i+3]])])
    data4_slice = module_copy.deepcopy(data4)
    data4_slice.reduce([select("time", times_hour[i])])
    g.plot(graph.data.list(data4_slice.data_center_list(strip_zero = True),
                           x = 1, y = 2,
                           title = "60-80"),
           styles = [graph.style.line(lineattrs = [color_list[i+4]])])
    data5_slice = module_copy.deepcopy(data5)
    data5_slice.reduce([select("time", times_hour[i])])
    g.plot(graph.data.list(data5_slice.data_center_list(strip_zero = True),
                           x = 1, y = 2,
                           title = "80-100"),
           styles = [graph.style.line(lineattrs = [color_list[i+5]])])


g.writePDFfile("out/testcase_withcoag/aero_comp_bc1_t6.pdf")
