#!/usr/bin/env python
# Copyright (C) 2007, 2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
sys.path.append("../tool")
from pmc_data_nc import *


subdir = "."
if len(sys.argv) > 1:
    subdir = sys.argv[1]

data = pmc_var(NetCDFFile("out/%s/urban_plume_0001.nc" % subdir),
	       "aero", [])

data.reduce([select("unit", "mass_den")])
data.write_summary(sys.stdout)
print
print "All values are in kg/m3"
print
for i_time in range(size(data.data,0)):
    print "time = ", data.dims[0].grid_centers[i_time], "seconds"
    print "dry_radius(m)",
    for i_spec in range(size(data.data,1)):
        print data.dims[1].grid_centers[i_spec],
    print
    for i_radius in range(size(data.data,2)):
        print data.dims[2].grid_centers[i_radius], 
        for i_spec in range(size(data.data,1)):
            print data.data[i_time,i_spec,i_radius],
        print
