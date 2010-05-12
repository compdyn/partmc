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
import config_filelist

data_base_dir = "data"
data_type = "diam_scrit_num"

x_axis = partmc.log_grid(min = config.diameter_axis_min,
                         max = config.diameter_axis_max,
                         n_bin = config.num_diameter_bins)
y_axis = partmc.log_grid(min = config.scrit_axis_min,
                         max = config.scrit_axis_max,
                         n_bin = config.num_scrit_bins)

def process_data(in_filename_list, out_filename):
    total_value = None
    for in_filename in in_filename_list:
        ncf = scipy.io.netcdf.netcdf_file(in_filename, 'r')
        particles = partmc.aero_particle_array_t(ncf)
        env_state = partmc.env_state_t(ncf)
        ncf.close()

        dry_diameters = particles.dry_diameters() * 1e6 # m to um
        scrit = (particles.critical_rel_humids(env_state) - 1) * 100 # in %

        value = partmc.histogram_2d(dry_diameters, scrit, x_axis, y_axis,
                                    weights = 1 / particles.comp_vols,
                                    only_positive=False)
        # do not account for y_axis in % as it is a log-scale
        value /= 1e6 # m^{-3} to cm^{-3}

        if total_value is None:
            total_value = value
        else:
            total_value += value
    total_value /= len(in_filename_list)
    np.savetxt(out_filename, total_value)
    mask = np.ma.make_mask(total_value <= 0.0)
    masked_total_value = np.ma.array(total_value, mask=mask)
    return(masked_total_value.min(), masked_total_value.max())

if __name__ == "__main__":
    global_min = None
    global_max = None
    for run in config_filelist.runs:
        data_dir = os.path.join(data_base_dir, run["name"])
        if not os.path.isdir(data_dir):
            os.mkdir(data_dir)
        for loop in run["loops"]:
            for index in loop["indices"]:
                data_name = "%s_%04d_%08d" % (data_type, loop["num"], index["num"])
                print run["name"] + " " + data_name
                in_filename_list = [proc["filename"] for proc in index["procs"]]
                out_filename = os.path.join(data_dir, data_name + ".txt")
                (value_min, value_max) = process_data(in_filename_list, out_filename)
                if global_min is None:
                    global_min = value_min
                else:
                    global_min = min(global_min, value_min)
                if global_max is None:
                    global_max = value_max
                else:
                    global_max = max(global_max, value_max)
                print "min = %g, max = %g" % (global_min, global_max)
