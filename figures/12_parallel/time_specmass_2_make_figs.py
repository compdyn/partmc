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
data_type = "time_specmass"

value_min = 1e-10
value_max = 1e-3

def make_plot(value, out_filename, aero_data):
    (figure, axes, cbar_axes) = config_matplotlib.make_fig(figure_width=6,
                                                           left_margin=0.8,
                                                           right_margin=2.4)

    axes.grid(True)
    #axes.grid(True, which = 'minor')
    axes.minorticks_on()

    xaxis = axes.get_xaxis()
    yaxis = axes.get_yaxis()
    xaxis.labelpad = 8
    yaxis.labelpad = 8

    axes.set_xlabel(r"elapsed time $t\ /\ \rm hr$")
    axes.set_ylabel(r"mass conc. $m\ /\ (\rm \mu g \ cm^{-3})$")

    handle_list = []
    label_list = []
    for i in range(1, np.size(value, 1)):
        h = plt.semilogy((value[:,0] - value[0,0]) / 3600, value[:,i])
        handle_list.append(h)
        name = aero_data.names[i - 1]
        label_list.append(partmc.aero_data_t.species_tex_names[name])
        plt.hold(True)
    axes.set_xbound(config.time_axis_min, config.time_axis_max)
    axes.set_ybound(value_min, value_max)

    plt.figlegend(handle_list, label_list, loc='center right', ncol=2)
    
    figure.savefig(out_filename)

if __name__ == "__main__":
    run = config_filelist.runs[0]
    loop = run["loops"][0]
    index = loop["indices"][0]
    proc = index["procs"][0]
    ncf = scipy.io.netcdf.netcdf_file(proc["filename"], 'r')
    aero_data = partmc.aero_data_t(ncf)
    ncf.close()
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
            make_plot(value, fig_filename, aero_data)
            plt.close('all')
