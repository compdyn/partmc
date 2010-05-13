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
data_type = "time_env"

value_temp_min = 285
value_temp_max = 305
value_rh_min = 50
value_rh_max = 100
value_pres_min = 90
value_pres_max = 110
value_height_min = 0
value_height_max = 500

def make_plot_temp(value, out_filename):
    (figure, axes, cbar_axes) = config_matplotlib.make_fig()

    axes.grid(True)
    #axes.grid(True, which = 'minor')
    axes.minorticks_on()

    xaxis = axes.get_xaxis()
    yaxis = axes.get_yaxis()
    xaxis.labelpad = 8
    yaxis.labelpad = 8

    axes.set_xlabel(r"elapsed time $t\ /\ \rm hr$")
    axes.set_ylabel(r"temperature $T\ /\ \rm K$")

    plt.plot((value[:,0] - value[0,0]) / 3600, value[:,1])
    axes.set_xbound(config.time_axis_min, config.time_axis_max)
    axes.set_ybound(value_temp_min, value_temp_max)

    figure.savefig(out_filename)

def make_plot_rh(value, out_filename):
    (figure, axes, cbar_axes) = config_matplotlib.make_fig()

    axes.grid(True)
    #axes.grid(True, which = 'minor')
    axes.minorticks_on()

    xaxis = axes.get_xaxis()
    yaxis = axes.get_yaxis()
    xaxis.labelpad = 8
    yaxis.labelpad = 8

    axes.set_xlabel(r"elapsed time $t\ /\ \rm hr$")
    axes.set_ylabel(r"relative humidity $\mbox{RH}\ /\ \rm \%$")

    plt.plot((value[:,0] - value[0,0]) / 3600, value[:,2] * 100)
    axes.set_xbound(config.time_axis_min, config.time_axis_max)
    axes.set_ybound(value_rh_min, value_rh_max)

    figure.savefig(out_filename)

def make_plot_pres(value, out_filename):
    (figure, axes, cbar_axes) = config_matplotlib.make_fig()

    axes.grid(True)
    #axes.grid(True, which = 'minor')
    axes.minorticks_on()

    xaxis = axes.get_xaxis()
    yaxis = axes.get_yaxis()
    xaxis.labelpad = 8
    yaxis.labelpad = 8

    axes.set_xlabel(r"elapsed time $t\ /\ \rm hr$")
    axes.set_ylabel(r"pressure $p\ /\ \rm kPa$")

    plt.plot((value[:,0] - value[0,0]) / 3600, value[:,3] / 1e3)
    axes.set_xbound(config.time_axis_min, config.time_axis_max )
    axes.set_ybound(value_pres_min, value_pres_max)

    figure.savefig(out_filename)

def make_plot_height(value, out_filename):
    (figure, axes, cbar_axes) = config_matplotlib.make_fig()

    axes.grid(True)
    #axes.grid(True, which = 'minor')
    axes.minorticks_on()

    xaxis = axes.get_xaxis()
    yaxis = axes.get_yaxis()
    xaxis.labelpad = 8
    yaxis.labelpad = 8

    axes.set_xlabel(r"elapsed time $t\ /\ \rm hr$")
    axes.set_ylabel(r"mixing layer height $H\ /\ \rm m$")

    plt.plot((value[:,0] - value[0,0]) / 3600, value[:,4])
    axes.set_xbound(config.time_axis_min, config.time_axis_max)
    axes.set_ybound(value_height_min, value_height_max)

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
            fig_filename_prefix = os.path.join(fig_dir, data_name + "_temp.pdf")
            make_plot_temp(value, fig_filename_prefix)
            fig_filename_prefix = os.path.join(fig_dir, data_name + "_rh.pdf")
            make_plot_rh(value, fig_filename_prefix)
            fig_filename_prefix = os.path.join(fig_dir, data_name + "_pres.pdf")
            make_plot_pres(value, fig_filename_prefix)
            fig_filename_prefix = os.path.join(fig_dir, data_name + "_height.pdf")
            make_plot_height(value, fig_filename_prefix)
            plt.close('all')
