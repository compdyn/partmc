#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from fig_helper import *
from numpy import *

out_file = "out/aging_aero_time_bulk.txt"

time_filename_list = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
n_time = len(time_filename_list)

data = zeros([n_time, 6])

for [i_time, [time, filename, key]] in enumerate(time_filename_list):
    print filename
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()

    num_den = 1.0 / array(particles.comp_vol)
    mass = particles.mass()
    dry_mass = particles.mass(exclude = ["H2O"])
    area = particles.surface_area()
    dry_area = particles.dry_surface_area()

    data[i_time, 0] = time

    data[i_time, 1] = num_den.sum()
    data[i_time, 2] = (num_den * mass).sum()
    data[i_time, 3] = (num_den * dry_mass).sum()
    data[i_time, 4] = (num_den * area).sum()
    data[i_time, 5] = (num_den * dry_area).sum()

savetxt(out_file, data)
