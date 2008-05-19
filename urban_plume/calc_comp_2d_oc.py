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

times_hour = [24]

netcdf_var = "comp_oc"
netcdf_dim = "composition_oc"
y_axis_label = r"$f_{{\rm BC},{\rm OC}}$"
filename = "figs/comp_2d_oc.pdf"

min_val = 0.0
max_val = 2.0

v_space = 0.5
h_space = 0.5

graph_width = 6

data = pmc_var(NetCDFFile("out/urban_plume_with_coag_0001.nc"),
	       netcdf_var,
	       [])
#data.write_summary(sys.stdout)

data.reduce([select("unit", "num_den"),
		 sum("aero_species")])
data.scale_dim(netcdf_dim, 100)
data.scale_dim("dry_radius", 2e6)
data.scale_dim("time", 1.0/3600)

data_slice = module_copy.deepcopy(data)
data_slice.reduce([select("time", times_hour[0])])

reducer1 = sum("dry_radius")
reducers = [reducer1]
data_slice.reduce(reducers)

data_help1 =  module_copy.deepcopy(data_slice)
data_help2 =  module_copy.deepcopy(data_slice)
data_help3 =  module_copy.deepcopy(data_slice)

reducer2 = sum("composition_oc", above = 90, below = 90)
reducers = [reducer2]
data_slice.reduce(reducers)
print data_slice.data

reducer3 = sum("composition_oc", above = 34, below = 34)
reducers = [reducer3]
data_help1.reduce(reducers)
print data_help1.data

reducer4 = sum("composition_oc", above = 0, below = 0)
reducers = [reducer4]
data_help2.reduce(reducers)
print data_help2.data

reducer5 = sum("composition_oc", above = 0, below = 100)
reducers = [reducer5]
data_help3.reduce(reducers)
print data_help3.data
print (data_help3.data - data_slice.data - data_help1.data - data_help2.data)*100/data_help3.data
