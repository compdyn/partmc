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
data_type = "time_mass"

value_min = 0
value_max = 0.0002

def make_plot(value, out_filename):
    (figure, axes, cbar_axes) = config_matplotlib.make_fig(left_margin=0.9)

    axes.grid(True)
    #axes.grid(True, which = 'minor')
    axes.minorticks_on()

    xaxis = axes.get_xaxis()
    yaxis = axes.get_yaxis()
    xaxis.labelpad = 8
    yaxis.labelpad = 8

    axes.set_xlabel(r"elapsed time $t\ /\ \rm hr$")
    axes.set_ylabel(r"mass conc. $m\ /\ (\rm \mu g \ cm^{-3})$")

    plt.plot((value[:,0] - value[0,0]) / 3600, value[:,1])
    axes.set_xbound(config.time_axis_min, config.time_axis_max)
    axes.set_ybound(value_min, value_max)
    
    figure.savefig(out_filename)

if __name__ == "__main__":
    for run in config_filelist.runs:
        data_dir = os.path.join(data_base_dir, run["name"])
        fig_dir = os.path.join(fig_base_dir, run["name"])
        if not os.path.isdir(fig_dir):
            os.mkdir(fig_dir)
        for loop in run["loops"]:
            data_name = "%s_%04d" % (data_type, loop["num"])
            print run["name"] + " " + data_name
            data_filename = os.path.join(data_dir, data_name + ".txt")
            value = np.loadtxt(data_filename)
            fig_filename = os.path.join(fig_dir, data_name + ".pdf")
            make_plot(value, fig_filename)
            plt.close('all')
