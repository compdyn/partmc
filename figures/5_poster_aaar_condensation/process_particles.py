#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math, re
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../../tool")
from pmc_data_nc import *
from pmc_pyx import *
import numpy
from config import *

data_prefix = "data/particles"

for [i_run, netcdf_pattern] in netcdf_indexed_patterns:
    time_filename_list = get_time_filename_list(netcdf_dir, netcdf_pattern)
    netcdf_filename = time_filename_list[0][1]
    ncf = NetCDFFile(netcdf_filename)
    particles = aero_particle_array_t(ncf)
    ncf.close()

    diameter = particles.diameter()

    init_grid = pmc_log_axis(min = 1e-8, max = 1e-6, n_bin = 20)
    ids = []
    for i in range(init_grid.n_bin):
        d = init_grid.center(i)
        i_p = (abs(diameter - d)).argmin()
        ids.append(particles.id[i_p])
    unique_ids = []
    for id in ids:
        if id not in unique_ids:
            unique_ids.append(id)
    ids = unique_ids
    print ids

    data = zeros([len(time_filename_list), len(ids) + 1])
    for (i_time, [time, filename, key]) in enumerate(time_filename_list):
        print filename
        ncf = NetCDFFile(filename)
        particles = aero_particle_array_t(ncf)
        env_state = env_state_t(ncf)
        ncf.close()

        data[i_time, 0] = env_state.elapsed_time
        diameter = particles.diameter()
        for (i_id, id) in enumerate(ids):
            i_p = (abs(particles.id - id)).argmin()
            if particles.id[i_p] != id:
                raise Exception("ID mismatch")
            data[i_time, i_id + 1] = diameter[i_p]

    data_filename = "%s_%d.txt" % (data_prefix, i_run)
    savetxt(data_filename, data)
