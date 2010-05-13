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
data_type = "diam_num"

value_min = 1
value_max = 35000

x_axis = partmc.log_grid(min = config.diameter_axis_min,
                         max = config.diameter_axis_max,
                         n_bin = config.num_diameter_bins)

def make_plot(value, out_filename):
    (figure, axes, cbar_axes) = config_matplotlib.make_fig()

    axes.grid(True)
    axes.grid(True, which = 'minor')
    axes.minorticks_on()
    axes.set_xscale('log')
    axes.set_yscale('log')

    xaxis = axes.get_xaxis()
    yaxis = axes.get_yaxis()
    xaxis.labelpad = 8
    yaxis.labelpad = 8

    #yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
    #yaxis.set_minor_locator(matplotlib.ticker.MaxNLocator(8))

    axes.set_xlabel(r"dry diameter $D_{\rm dry}\ /\ \rm\mu m$")
    axes.set_ylabel(r"number conc. $n\ /\ \rm cm^{-3}$")

    plt.loglog(x_axis.centers(), value)
    axes.set_xbound(x_axis.min, x_axis.max)
    axes.set_ybound(value_min, value_max)
    
    figure.savefig(out_filename)

if __name__ == "__main__":
    for run in config_filelist.runs:
        data_dir = os.path.join(data_base_dir, run["name"])
        fig_dir = os.path.join(fig_base_dir, run["name"])
        if not os.path.isdir(fig_dir):
            os.mkdir(fig_dir)
        for loop in run["loops"]:
            for index in loop["indices"]:
                data_name = "%s_%04d_%08d" % (data_type, loop["num"], index["num"])
                print run["name"] + " " + data_name
                data_filename = os.path.join(data_dir, data_name + ".txt")
                value = np.loadtxt(data_filename)
                fig_filename = os.path.join(fig_dir, data_name + ".pdf")
                make_plot(value, fig_filename)
                plt.close('all')
