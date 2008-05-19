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

data_no = pmc_var(NetCDFFile("out/urban_plume_no_coag_0001.nc"),
	       "comp_bc",
	       [])
data_wc = pmc_var(NetCDFFile("out/urban_plume_with_coag_0001.nc"),
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

for i in range(len(composition_lower)):
    data_no_slice = module_copy.deepcopy(data_no)
    data_wc_slice = module_copy.deepcopy(data_wc)
        
    reducer1 = sum("composition_bc", above = composition_lower[i], below = composition_upper[i])
    reducers = [reducer1]
    data_no_slice.reduce(reducers)
    data_wc_slice.reduce(reducers)

    data_no_total1 = module_copy.deepcopy(data_no_slice)
    data_wc_total1 = module_copy.deepcopy(data_wc_slice)

    data_no_total2 = module_copy.deepcopy(data_no_slice)
    data_wc_total2 = module_copy.deepcopy(data_wc_slice)
    
    reducer2 = sum("dry_radius", below = 0.04)
    reducers = [reducer2]
    data_no_total1.reduce(reducers)
    data_wc_total1.reduce(reducers)

    print i, data_no_total1.data, data_wc_total1.data, (data_wc_total1.data - data_no_total1.data)*100/data_no_total1.data

    reducer3 = sum("dry_radius", above = 0.1)
    reducers = [reducer3]
    data_no_total2.reduce(reducers)
    data_wc_total2.reduce(reducers)

    print i, data_no_total2.data, data_wc_total2.data, (data_wc_total2.data - data_no_total2.data)*100/data_no_total2.data

