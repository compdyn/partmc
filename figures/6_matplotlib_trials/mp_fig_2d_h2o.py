#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import numpy as np
sys.path.append("../../tool")
from pmc_data_nc import *
from config import *

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

const = load_constants("../../src/constants.f90")

out_prefix = "figs/mp_2d_h2o"

def get_plot_data_bc(filename, value_min = None, value_max = None):
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()

    diameter = particles.dry_diameter() * 1e6
    comp_frac = particles.mass(include = ["BC"]) \
                / particles.mass(exclude = ["H2O"]) * 100
    # hack to avoid landing just around the integer boundaries
    comp_frac *= (1.0 + 1e-12)
    h2o = particles.mass(include = ["H2O"])

    x_axis = pmc_log_axis(min = diameter_axis_min, max = diameter_axis_max,
                          n_bin = num_diameter_bins * 2)
    y_axis = pmc_linear_axis(min = bc_axis_min, max = bc_axis_max,
                             n_bin = num_bc_bins * 2)

    value = multival_2d(diameter, comp_frac, h2o, x_axis, y_axis)
    if value_max == None:
        value_max = value.max()
    if value_min == None:
        maxed_value = where(value > 0.0, value, value_max)
        value_min = maxed_value.min()
    #if value_max > 0.0:
    #    value = (log(value) - log(value_min)) \
    #            / (log(value_max) - log(value_min))
    #value = value.clip(0.0, 1.0)

    return (value, x_axis.edges(), y_axis.edges(),
            env_state, value_min, value_max)

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

def make_2d_plot(in_filename, out_filename):
    (figure, axes, cbar_axes) = make_fig(colorbar = True, right_margin = 0.9)

    axes.grid(True)
    axes.grid(True, which = 'minor')
    axes.minorticks_on()
    axes.set_xscale('log')

    #axes.set_xticks([0, 6, 12, 18, 24])
    #axes.set_xticks([3, 9, 15, 21], minor = True)
    #axes.set_yticks([0.025, 0.075, 0.125, 0.175], minor = True)
    
    axes.set_xbound(diameter_axis_min, diameter_axis_max)
    axes.set_ybound(bc_axis_min, bc_axis_max)

    xaxis = axes.get_xaxis()
    yaxis = axes.get_yaxis()
    xaxis.labelpad = 8
    yaxis.labelpad = 8
    #xaxis.set_major_formatter(matplotlib.ticker.LogFormatter())
    yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
    yaxis.set_minor_locator(matplotlib.ticker.MaxNLocator(8))

    axes.set_xlabel(r"dry diameter $D\ (\rm\mu m)$")
    axes.set_ylabel(r"BC dry mass frac. $w_{{\rm BC},{\rm dry}}\ (\%)$")

    (value, x_edges, y_edges, env_state, value_min, value_max) = get_plot_data_bc(in_filename)

    axes.set_xbound(diameter_axis_min, diameter_axis_max)
    axes.set_ybound(bc_axis_min, bc_axis_max)
    
    p = axes.pcolor(x_edges, y_edges, value.transpose(), norm = matplotlib.colors.LogNorm(),
                    cmap=matplotlib.cm.jet, linewidths = 0.1)
    figure.colorbar(p, cax = cbar_axes, format = matplotlib.ticker.LogFormatterMathtext())
    cbar_axes.set_ylabel(r"water mass $(\rm kg)$")
    cbar_yaxis = cbar_axes.get_yaxis()

    figure.savefig(out_filename)

for [i_run, netcdf_pattern] in netcdf_indexed_patterns:
    out_filename = "%s_%d.pdf" % (out_prefix, i_run)
    print out_filename

    filename_list = get_filename_list(netcdf_dir, netcdf_pattern)
    in_filename = filename_list[0]
    make_2d_plot(in_filename, out_filename)
    
