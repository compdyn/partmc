#!/usr/bin/env python
# Copyright (C) 2007-2010 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import numpy as np
import scipy.io
sys.path.append("../../tool")
import partmc
import config

data_base_dir = "data"
data_type = "2d_bc"

x_axis = partmc.log_grid(min = config.diameter_axis_min,
                         max = config.diameter_axis_max,
                         n_bin = config.num_diameter_bins)
y_axis = partmc.linear_grid(min = config.bc_axis_min,
                            max = config.bc_axis_max,
                            n_bin = config.num_bc_bins)

for run in config.runs:
    data_dir = os.path.join(data_base_dir, run["name"])
    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)
    for loop in run["loops"]:
        for index in loop["indices"]:
            data_name = "%s_%04d_%08d" % (data_type, loop["num"], index["num"])
            print run["name"] + " " + data_name
            total_value = None
            for proc in index["procs"]:
                ncf = scipy.io.netcdf.netcdf_file(proc["filename"], 'r')
                particles = partmc.aero_particle_array_t(ncf)
                env_state = partmc.env_state_t(ncf)
                ncf.close()

                diameters = particles.dry_diameters() * 1e6
                comp_frac = particles.masses(include = ["BC"]) \
                            / particles.masses(exclude = ["H2O"]) * 100

                # hack to avoid landing just around the integer boundaries
                comp_frac *= (1.0 + 1e-12)

                value = partmc.histogram_2d(diameters, comp_frac, x_axis, y_axis,
                                            weights = 1 / particles.comp_vols)
                value *= 100
                value /= 1e6
                
                if total_value is None:
                    total_value = value
                else:
                    total_value += value
            total_value /= len(index["procs"])
            data_filename = os.path.join(data_dir, data_name + ".txt")
            np.savetxt(data_filename, total_value)
