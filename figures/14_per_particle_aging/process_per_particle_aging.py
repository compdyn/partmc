#!/usr/bin/env python

import os, sys
import scipy.io
import numpy as np

sys.path.append("/Users/nriemer/subversion/partmc/trunk/tool")
import partmc

class Struct(object): 
	def __init__(self): 
		pass

particle_set = {}

netcdf_dir = "/Users/nriemer/subversion/partmc/trunk/scenarios/2_urban_plume2/out"
netcdf_pattern = "urban_plume_nc_0001_(.*).nc"
time_filename_list = partmc.get_time_filename_list(netcdf_dir, netcdf_pattern)

for [time, filename, key] in time_filename_list:
	print time, filename, key
	ncf = scipy.io.netcdf.netcdf_file(filename, 'r')
	particles = partmc.aero_particle_array_t(ncf)
	env_state = partmc.env_state_t(ncf)
	ncf.close()

