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
data.scale_dim("dry_radius", 2e6)
data.scale_dim("time", 1.0/3600)

data.dim_by_name("n_orig_part").grid_centers = data.dim_by_name("n_orig_part").grid_centers - 1
data.dim_by_name("n_orig_part").grid_edges = data.dim_by_name("n_orig_part").grid_edges - 1

data_slice = module_copy.deepcopy(data)

data_slice.reduce([select("time", time_hour)])
reducer2 = sum("dry_radius")
reducers = [reducer2]
data_slice.reduce(reducers)

data_help = module_copy.deepcopy(data_slice)

reducer3 = sum("n_orig_part")
reducers = [reducer3]
data_slice.reduce(reducers)

reducer4 = sum("n_orig_part", above = 6) # means 5 coag events or more
reducers = [reducer4]
data_help.reduce(reducers)

print data_slice.data, data_help.data, data_help.data/data_slice.data
