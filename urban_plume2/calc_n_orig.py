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

# fix data to center on integers
data.dim_by_name("n_orig_part").grid_centers -= 0.5
data.dim_by_name("n_orig_part").grid_edges -= 0.5

# shift from "num constituent particles" to "num coag events"
data.dim_by_name("n_orig_part").grid_centers -= 1
data.dim_by_name("n_orig_part").grid_edges -= 1

data_slice = module_copy.deepcopy(data)
data_slice.reduce([select("time", time_hour),
                   sum("dry_radius")])
data_help = module_copy.deepcopy(data_slice)
data_slice.reduce([sum("n_orig_part")])

print "total =", data_slice.data

print "%20s %20s %20s" % ("num_events", "num_den", "percentage")
for i in range(30):
    data_help2 = module_copy.deepcopy(data_help)
    data_help2.reduce([sum("n_orig_part", above = i)]) # >= 5 coag events
    print "%20i %20f %20.2f%%" % (i, data_help2.data,
                                  data_help2.data / data_slice.data * 100)
