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
data_type = "time_env"

def process_data(in_filename_list):
    total_value = None
    for in_filename in in_filename_list:
        ncf = scipy.io.netcdf.netcdf_file(in_filename, 'r')
        env_state = partmc.env_state_t(ncf)
        ncf.close()

        time_since_midnight = env_state.start_time_of_day \
            + env_state.elapsed_time

        value = np.array([env_state.temperature,
                          env_state.relative_humidity,
                          env_state.pressure,
                          env_state.height,
                          ])

        if total_value is None:
            total_value = np.array(value)
        else:
            total_value += value
    total_value /= len(in_filename_list)
    return np.concatenate((np.array([time_since_midnight]), total_value))

if __name__ == "__main__":
    for run in config_filelist.runs:
        data_dir = os.path.join(data_base_dir, run["name"])
        if not os.path.isdir(data_dir):
            os.mkdir(data_dir)
        for loop in run["loops"]:
            data_name = "%s_%04d" % (data_type, loop["num"])
            print run["name"] + " " + data_name
            out_filename = os.path.join(data_dir, data_name + ".txt")
            for (i_index, index) in enumerate(loop["indices"]):
                in_filename_list = [proc["filename"] for proc in index["procs"]]
                index_value = process_data(in_filename_list)
                if i_index == 0:
                    value = np.zeros((len(loop["indices"]), len(index_value)))
                value[i_index,:] = index_value
            np.savetxt(out_filename, value)
