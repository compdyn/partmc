#!/usr/bin/env python
# Copyright (C) 2007-2010 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import numpy as np
import scipy.io
sys.path.append("../../tool")
import partmc
import config

matplotlib.rc('text', usetex = True)
matplotlib.rc('xtick.major', pad = 8)
matplotlib.rc('ytick.major', pad = 8)
matplotlib.rc('xtick', labelsize = 10)
matplotlib.rc('legend', fontsize = 10, borderpad = 0.7, borderaxespad = 1)
matplotlib.rc('font', size = 10, family = "serif",
              serif = ["Computer Modern Roman"])
matplotlib.rc('lines', linewidth = 0.5)
matplotlib.rc('patch', linewidth = 0.5)
matplotlib.rc('axes', linewidth = 0.5)

fig_base_dir = "figs"
data_base_dir = "data"
data_type = "2d_bc"

x_axis = partmc.log_grid(min = config.diameter_axis_min,
                         max = config.diameter_axis_max,
                         n_bin = config.num_diameter_bins)
y_axis = partmc.linear_grid(min = config.bc_axis_min,
                            max = config.bc_axis_max,
                            n_bin = config.num_bc_bins)

def make_fig(figure_width = 4,
             figure_height = None,
             axis_ratio = (1 + math.sqrt(5)) / 2, # golden ratio
             left_margin = 0.6,
             right_margin = 0.2,
             bottom_margin = 0.5,
             top_margin = 0.2,
             colorbar = False,
             colorbar_width = 0.15,
             colorbar_height_fraction = 0.8,
             colorbar_offset = 0.2):
    axis_width = figure_width - left_margin - right_margin
    axis_height = axis_width / axis_ratio
    figure_height = bottom_margin + axis_height + top_margin
    left_margin_fraction = left_margin / figure_width
    bottom_margin_fraction = bottom_margin / figure_height
    axis_width_fraction = axis_width / figure_width
    axis_height_fraction = axis_height / figure_height
    figure = plt.figure()
    figure.set_figwidth(figure_width)
    figure.set_figheight(figure_height)
    axes = figure.add_axes([left_margin_fraction,
                            bottom_margin_fraction,
                            axis_width_fraction,
                            axis_height_fraction])
    if colorbar:
        cb_left_fraction = (left_margin + axis_width + colorbar_offset) / figure_width
        cb_bottom_fraction = (bottom_margin + axis_height * (1.0 - colorbar_height_fraction) / 2.0) / figure_height
        cb_width_fraction = colorbar_width / figure_width
        cb_height_fraction = axis_height * colorbar_height_fraction / figure_height
        colorbar_axes = figure.add_axes([cb_left_fraction,
                                         cb_bottom_fraction,
                                         cb_width_fraction,
                                         cb_height_fraction])
    else:
        colorbar_axes = None
    return (figure, axes, colorbar_axes)

def make_2d_plot(value, out_filename):
    (figure, axes, cbar_axes) = make_fig(colorbar = True, right_margin = 0.9)

    axes.grid(True)
    axes.grid(True, which = 'minor')
    axes.minorticks_on()
    axes.set_xscale('log')

    axes.set_xbound(config.diameter_axis_min, config.diameter_axis_max)
    axes.set_ybound(config.bc_axis_min, config.bc_axis_max)

    xaxis = axes.get_xaxis()
    yaxis = axes.get_yaxis()
    xaxis.labelpad = 8
    yaxis.labelpad = 8

    yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
    yaxis.set_minor_locator(matplotlib.ticker.MaxNLocator(8))

    axes.set_xlabel(r"dry diameter $D\ (\rm\mu m)$")
    axes.set_ylabel(r"BC dry mass frac. $w_{{\rm BC},{\rm dry}}\ (\%)$")

    axes.set_xbound(config.diameter_axis_min, config.diameter_axis_max)
    axes.set_ybound(config.bc_axis_min, config.bc_axis_max)
    
    p = axes.pcolor(x_axis.edges(), y_axis.edges(), value.transpose(),
                    norm = matplotlib.colors.LogNorm(),
                    cmap=matplotlib.cm.jet, linewidths = 0.1)
    figure.colorbar(p, cax = cbar_axes, format = matplotlib.ticker.LogFormatterMathtext())
    cbar_axes.set_ylabel(r"number conc. $(\rm cm^{-3})$")
    figure.savefig(out_filename)

for run in config.runs:
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
            mask = np.ma.make_mask(value <= 0.0)
            value = np.ma.array(value, mask=mask)
            fig_filename = os.path.join(fig_dir, data_name + ".pdf")
            make_2d_plot(value, fig_filename)
