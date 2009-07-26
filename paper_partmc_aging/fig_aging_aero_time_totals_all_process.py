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

out_file = "out/aging_aero_time_totals_all.txt"

const = load_constants("../src/constants.f90")

time_filename_list = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
n_time = len(time_filename_list)

data = zeros([n_time, 10])

for [i_time, [time, filename, key]] in enumerate(time_filename_list):
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()

    num_den = 1.0 / array(particles.comp_vol)
    mass = particles.mass()
    area = particles.surface_area()
    soot_mass = particles.mass(include = ["BC"])
    critical_ss = particles.kappa_rh(env_state, const) - 1.0
    is_aged = (critical_ss < level_mid_value)
    is_fresh = logical_not(is_aged)

    num_den_bc = where(soot_mass > 0, num_den, 0)
    num_den_bc_aged = where(is_aged, num_den_bc, 0)
    num_den_bc_fresh = where(is_fresh, num_den_bc, 0)

    data[i_time, 0] = time

    data[i_time, 1] = num_den_bc.sum()
    data[i_time, 2] = num_den_bc_aged.sum()
    data[i_time, 3] = num_den_bc_fresh.sum()
    
    data[i_time, 4] = (num_den_bc * mass).sum()
    data[i_time, 5] = (num_den_bc_aged * mass).sum()
    data[i_time, 6] = (num_den_bc_fresh * mass).sum()
    
    data[i_time, 7] = (num_den_bc * area).sum()
    data[i_time, 8] = (num_den_bc_aged * area).sum()
    data[i_time, 9] = (num_den_bc_fresh * area).sum()

    print i_time, data[i_time, 0], data[i_time, 1], data[i_time, 4], data[i_time, 7]

savetxt(out_file, data)
