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

value_min = 80
value_max = 1e6

def make_plot(error_data, out_filename):
    (figure, axes, cbar_axes) = config_matplotlib.make_fig()

    axes.grid(True)
    #axes.grid(True, which = 'minor')
    axes.minorticks_on()
    axes.set_xscale('log')
    axes.set_yscale('log')

    xaxis = axes.get_xaxis()
    yaxis = axes.get_yaxis()
    xaxis.labelpad = 8
    yaxis.labelpad = 8

    axes.set_xlabel(r"number of cores")
    axes.set_ylabel(r"$L_2$ error")

    #axes.set_xbound(x_axis.min, x_axis.max)
    #axes.set_ybound(y_axis.min, y_axis.max)

    plt.loglog(error_data[:,0], error_data[:,1])
    figure.savefig(out_filename)

if __name__ == "__main__":
    for run_group in config_filelist.runs_by_base:
        fig_dir = os.path.join(fig_base_dir, run_group["name"])
        if not os.path.isdir(fig_dir):
            os.mkdir(fig_dir)
        error_data_dir = os.path.join(data_base_dir, run_group["name"])
        if not os.path.isdir(error_data_dir):
            os.mkdir(error_data_dir)
        first_run = run_group["run_items"][0]["run"]
        for (i_loop, first_loop) in enumerate(first_run["loops"]):
            for (i_index, first_index) in enumerate(first_loop["indices"]):
                data_name = "%s_%04d_%08d" % (data_type, first_loop["num"], first_index["num"])
                print run_group["name"] + " " + data_name
                error_data = np.zeros((len(run_group["run_items"]), 2))
                for (i_run_item, run_item) in enumerate(run_group["run_items"]):
                    run = run_item["run"]
                    data_dir = os.path.join(data_base_dir, run["name"])
                    loop = run["loops"][i_loop]
                    index = loop["indices"][i_index]
                    error_filename = os.path.join(data_dir, data_name + "_error.txt")
                    error = np.loadtxt(error_filename)
                    error_data[i_run_item, 0] = run_item["size"]
                    error_data[i_run_item, 1] = float(error)
                error_data_filename = os.path.join(error_data_dir, data_name + "_error.txt")
                np.savetxt(error_data_filename, error_data)
                fig_filename = os.path.join(fig_dir, data_name + "_error.pdf")
                make_plot(error_data, fig_filename)
                plt.close('all')
