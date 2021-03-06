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
data_type = "time_num"

def process_data(in_filename_list):
    total_value = None
    for in_filename in in_filename_list:
        ncf = scipy.io.netcdf.netcdf_file(in_filename, 'r')
        particles = partmc.aero_particle_array_t(ncf)
        env_state = partmc.env_state_t(ncf)
        ncf.close()

        time_since_midnight = env_state.start_time_of_day \
            + env_state.elapsed_time

        value = (1 / particles.comp_vols).sum()
        value /= 1e6 # m^{-3} to cm^{-3}

        if total_value is None:
            total_value = value
        else:
            total_value += value
    total_value /= len(in_filename_list)
    return [time_since_midnight, total_value]

if __name__ == "__main__":
    global_min = None
    global_max = None
    for run in config_filelist.runs:
        data_dir = os.path.join(data_base_dir, run["name"])
        if not os.path.isdir(data_dir):
            os.mkdir(data_dir)
        for loop in run["loops"]:
            data_name = "%s_%04d" % (data_type, loop["num"])
            print run["name"] + " " + data_name
            out_filename = os.path.join(data_dir, data_name + ".txt")
            value = np.zeros((len(loop["indices"]), 2))
            for (i_index, index) in enumerate(loop["indices"]):
                in_filename_list = [proc["filename"] for proc in index["procs"]]
                value[i_index,:] = process_data(in_filename_list)
            np.savetxt(out_filename, value)
            if global_min is None:
                global_min = value[:,1].min()
            else:
                global_min = min(global_min, value[:,1].min())
            if global_max is None:
                global_max = value[:,1].max()
            else:
                global_max = max(global_max, value[:,1].max())
            print "min = %g, max = %g" % (global_min, global_max)
