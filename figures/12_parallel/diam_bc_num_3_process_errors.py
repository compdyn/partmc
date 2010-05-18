#!/usr/bin/env python
# Copyright (C) 2007-2010 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import config
import config_filelist
import config_matplotlib
import os, sys, math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
sys.path.append("../../tool")
import partmc

fig_base_dir = "figs"
data_base_dir = "data"
data_type = "diam_bc_num"

x_axis = partmc.log_grid(min = config.diameter_axis_min,
                         max = config.diameter_axis_max,
                         n_bin = config.num_diameter_bins)
y_axis = partmc.linear_grid(min = config.bc_axis_min,
                            max = config.bc_axis_max,
                            n_bin = config.num_bc_bins)

def compute_error(value, true_value, out_filename):
    error = (value - true_value) * x_axis.grid_size(0) * y_axis.grid_size(0)
    error_norm = np.sqrt((error**2).sum())
    np.savetxt(out_filename, np.array([error_norm]))

if __name__ == "__main__":
    true_run = config_filelist.true_run
    for run in config_filelist.runs:
        data_dir = os.path.join(data_base_dir, run["name"])
        true_data_dir = os.path.join(data_base_dir, true_run["name"])
        fig_dir = os.path.join(fig_base_dir, run["name"])
        if not os.path.isdir(fig_dir):
            os.mkdir(fig_dir)
        for (loop, true_loop) in zip(run["loops"], true_run["loops"]):
            for (index, true_index) in zip(loop["indices"], true_loop["indices"]):
                data_name = "%s_%04d_%08d" % (data_type, loop["num"], index["num"])
                true_data_name = "%s_%04d_%08d" % (data_type, true_loop["num"], true_index["num"])
                print run["name"] + " " + data_name
                data_filename = os.path.join(data_dir, data_name + ".txt")
                true_data_filename = os.path.join(true_data_dir, true_data_name + ".txt")
                value = np.loadtxt(data_filename)
                true_value = np.loadtxt(true_data_filename)
                out_filename = os.path.join(data_dir, data_name + "_error.txt")
                compute_error(value, true_value, out_filename)
