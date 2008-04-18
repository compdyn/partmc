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

times_hour = [24]

data1 = pmc_var(NetCDFFile("out/testcase_nocoag/urban_plume_state_0001.nc"),
	       "comp_bc",
	       [])
data1.write_summary(sys.stdout)
data1.scale_dim("composition", 100)
data1.scale_dim("radius", 1e6)
data1.scale_dim("time", 1.0/3600)

data1.reduce([select("unit", "vol_den"),
		 sum("aero_species"),
                 sum("composition", above = 0, below = 2)])

data2 = pmc_var(NetCDFFile("out/testcase_nocoag/urban_plume_state_0001.nc"),
               "comp_bc",
               [])
data2.write_summary(sys.stdout)
data2.scale_dim("composition", 100)
data2.scale_dim("radius", 1e6)
data2.scale_dim("time", 1.0/3600)

data2.reduce([select("unit", "vol_den"),
                 sum("aero_species"),
                 sum("composition", above = 2, below = 10)])

data3 = pmc_var(NetCDFFile("out/testcase_nocoag/urban_plume_state_0001.nc"),
               "comp_bc",
               [])
data3.write_summary(sys.stdout)
data3.scale_dim("composition", 100)
data3.scale_dim("radius", 1e6)
data3.scale_dim("time", 1.0/3600)

data3.reduce([select("unit", "vol_den"),
                 sum("aero_species"),
                 sum("composition", above = 10, below = 20)])

data4 = pmc_var(NetCDFFile("out/testcase_nocoag/urban_plume_state_0001.nc"),
               "comp_bc",
               [])
data4.write_summary(sys.stdout)
data4.scale_dim("composition", 100)
data4.scale_dim("radius", 1e6)
data4.scale_dim("time", 1.0/3600)

data4.reduce([select("unit", "vol_den"),
                 sum("aero_species"),
                 sum("composition", above = 20, below = 40)])

data5 = pmc_var(NetCDFFile("out/testcase_nocoag/urban_plume_state_0001.nc"),
               "comp_bc",
               [])
data5.write_summary(sys.stdout)
data5.scale_dim("composition", 100)
data5.scale_dim("radius", 1e6)
data5.scale_dim("time", 1.0/3600)

data5.reduce([select("unit", "vol_den"),
                 sum("aero_species"),
                 sum("composition", above = 80, below = 90)])

data1a = pmc_var(NetCDFFile("out/testcase_withcoag/urban_plume_state_0001.nc"),
	       "comp_bc",
	       [])
data1a.write_summary(sys.stdout)
data1a.scale_dim("composition", 100)
data1a.scale_dim("radius", 1e6)
data1a.scale_dim("time", 1.0/3600)

data1a.reduce([select("unit", "vol_den"),
		 sum("aero_species"),
                 sum("composition", above = 0, below = 2)])

data2a = pmc_var(NetCDFFile("out/testcase_withcoag/urban_plume_state_0001.nc"),
               "comp_bc",
               [])
data2a.write_summary(sys.stdout)
data2a.scale_dim("composition", 100)
data2a.scale_dim("radius", 1e6)
data2a.scale_dim("time", 1.0/3600)

data2a.reduce([select("unit", "vol_den"),
                 sum("aero_species"),
                 sum("composition", above = 2, below = 10)])

data3a = pmc_var(NetCDFFile("out/testcase_withcoag/urban_plume_state_0001.nc"),
               "comp_bc",
               [])
data3a.write_summary(sys.stdout)
data3a.scale_dim("composition", 100)
data3a.scale_dim("radius", 1e6)
data3a.scale_dim("time", 1.0/3600)

data3a.reduce([select("unit", "vol_den"),
                 sum("aero_species"),
                 sum("composition", above = 10, below = 20)])

data4a = pmc_var(NetCDFFile("out/testcase_withcoag/urban_plume_state_0001.nc"),
               "comp_bc",
               [])
data4a.write_summary(sys.stdout)
data4a.scale_dim("composition", 100)
data4a.scale_dim("radius", 1e6)
data4a.scale_dim("time", 1.0/3600)

data4a.reduce([select("unit", "vol_den"),
                 sum("aero_species"),
                 sum("composition", above = 20, below = 40)])

data5a = pmc_var(NetCDFFile("out/testcase_withcoag/urban_plume_state_0001.nc"),
               "comp_bc",
               [])
data5a.write_summary(sys.stdout)
data5a.scale_dim("composition", 100)
data5a.scale_dim("radius", 1e6)
data5a.scale_dim("time", 1.0/3600)

data5a.reduce([select("unit", "vol_den"),
                 sum("aero_species"),
                 sum("composition", above = 80, below = 90)])


g = graph.graphxy(
	width = 10,
	x = graph.axis.log(title = r'radius ($\mu$m)',
                           min = 0.01, max = 1,
			   painter = grid_painter),
	y = graph.axis.log(min = 1e-16, max = 1e-9,
			      title = r'volume density ($\rm m^{-3} \, m^{-3}$)',
			      painter = grid_painter),
        key = graph.key.key(pos = "tr"))

for i in range(len(times_hour)):
    data1_slice = module_copy.deepcopy(data1)
    data1_slice.reduce([select("time", times_hour[i])])
    g.plot(graph.data.list(data1_slice.data_center_list(strip_zero = True),
			   x = 1, y = 2, 
                           title = "0-2 no"),
	   styles = [graph.style.line(lineattrs = [color_list[i+1]])])
    data2_slice = module_copy.deepcopy(data2)
    data2_slice.reduce([select("time", times_hour[i])])
    g.plot(graph.data.list(data2_slice.data_center_list(strip_zero = True),
                           x = 1, y = 2,
                           title = "2-10 no"),
           styles = [graph.style.line(lineattrs = [color_list[i+2]])])
    data3_slice = module_copy.deepcopy(data3)
    data3_slice.reduce([select("time", times_hour[i])])
#    g.plot(graph.data.list(data3_slice.data_center_list(strip_zero = True),
#                           x = 1, y = 2,
#                           title = "10-20 no"),
#           styles = [graph.style.line(lineattrs = [color_list[i+3]])])
    data4_slice = module_copy.deepcopy(data4)
    data4_slice.reduce([select("time", times_hour[i])])
#    g.plot(graph.data.list(data4_slice.data_center_list(strip_zero = True),
#                           x = 1, y = 2,
#                           title = "20-40"),
#           styles = [graph.style.line(lineattrs = [color_list[i+4]])])
    data5_slice = module_copy.deepcopy(data5)
    data5_slice.reduce([select("time", times_hour[i])])
    g.plot(graph.data.list(data5_slice.data_center_list(strip_zero = True),
                           x = 1, y = 2,
                           title = "80-90 no"),
           styles = [graph.style.line(lineattrs = [color_list[i+5]])])
    data1a_slice = module_copy.deepcopy(data1a)
    data1a_slice.reduce([select("time", times_hour[i])])
    g.plot(graph.data.list(data1a_slice.data_center_list(strip_zero = True),
			   x = 1, y = 2, 
                           title = "0-2 with"),
	   styles = [graph.style.line(lineattrs = [color_list[i+1],style.linewidth.THick])])
    data2a_slice = module_copy.deepcopy(data2a)
    data2a_slice.reduce([select("time", times_hour[i])])
    g.plot(graph.data.list(data2a_slice.data_center_list(strip_zero = True),
                           x = 1, y = 2,
                           title = "2-10 with"),
           styles = [graph.style.line(lineattrs = [color_list[i+2],style.linewidth.THick])])
    data3a_slice = module_copy.deepcopy(data3)
    data3a_slice.reduce([select("time", times_hour[i])])
#    g.plot(graph.data.list(data3a_slice.data_center_list(strip_zero = True),
#                           x = 1, y = 2,
#                           title = "10-20 with"),
#           styles = [graph.style.line(lineattrs = [color_list[i+3],style.linewidth.THick])])
    data4a_slice = module_copy.deepcopy(data4a)
    data4a_slice.reduce([select("time", times_hour[i])])
#    g.plot(graph.data.list(data4a_slice.data_center_list(strip_zero = True),
#                           x = 1, y = 2,
#                           title = "20-40"),
#           styles = [graph.style.line(lineattrs = [color_list[i+4]])])
    data5a_slice = module_copy.deepcopy(data5a)
    data5a_slice.reduce([select("time", times_hour[i])])
    g.plot(graph.data.list(data5a_slice.data_center_list(strip_zero = True),
                           x = 1, y = 2,
                           title = "80-90 with"),
           styles = [graph.style.line(lineattrs = [color_list[i+5],style.linewidth.THick])])


g.writePDFfile("out/testcase_nocoag/aero_comp_bc1v_t24ww.pdf")
