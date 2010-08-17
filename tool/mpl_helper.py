#!/usr/bin/env python
# Copyright (C) 2007-2010 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import sys, os
import math
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

matplotlib.rc('text', usetex = True)
matplotlib.rc('xtick.major', pad = 8)
matplotlib.rc('ytick.major', pad = 8)
matplotlib.rc('xtick', labelsize = 10)
matplotlib.rc('legend', fontsize = 10, borderpad = 0.7, borderaxespad = 1)
matplotlib.rc('font', size = 10, family = "serif",
              serif = ["Computer Modern Roman"])
matplotlib.rc('lines', linewidth = 1.5)
matplotlib.rc('patch', linewidth = 0.5)
matplotlib.rc('axes', linewidth = 0.5)

def make_fig(figure_width=5,
             axis_ratio=(1 + math.sqrt(5)) / 2, # golden ratio
             left_margin=0.8,
             right_margin=0.2,
             bottom_margin=0.5,
             top_margin=0.2,
             colorbar=False,
             colorbar_width=0.15,
             colorbar_height_fraction=0.8,
             colorbar_offset=0.2,
             ):
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
        return (figure, axes, colorbar_axes)
    else:
        return (figure, axes)

def make_fig_array(n_vert=2,
                   n_horiz=2,
                   figure_width=5,
                   axis_ratio=(1 + math.sqrt(5)) / 2, # golden ratio
                   left_margin=0.8,
                   right_margin=0.2,
                   bottom_margin=0.5,
                   top_margin=0.2,
                   horiz_sep=0.4,
                   vert_sep=0.4,
                   colorbar=False,
                   colorbar_width=0.15,
                   colorbar_height_fraction=0.8,
                   colorbar_offset=0.2,
                   ):
    """
    Numbering convention:
    axes_array[i][j] has (i,j):

          (2,0) (2,1)
          (1,0) (1,1)
          (0,0) (0,1)
    """
    axis_width = (figure_width - left_margin - right_margin - (n_horiz - 1) * horiz_sep) / float(n_horiz)
    axis_height = axis_width / axis_ratio
    figure_height = bottom_margin + axis_height * n_vert + vert_sep * (n_vert - 1) + top_margin
    figure = plt.figure()
    figure.set_figwidth(figure_width)
    figure.set_figheight(figure_height)
    axes_array = []
    for i_vert in range(n_vert):
        axes_array.append([])
        for i_horiz in range(n_horiz):
            x_left = left_margin + i_horiz * (axis_width + horiz_sep)
            y_bottom = bottom_margin + i_vert * (axis_height + vert_sep)
            kwargs = {}
            if i_horiz > 0:
                kwargs["sharey"] = last_y_axes
            if i_vert > 0:
                kwargs["sharex"] = last_x_axes
            new_axes = figure.add_axes([x_left / figure_width,
                                        y_bottom / figure_height,
                                        axis_width / figure_width,
                                        axis_height / figure_height],
                                       **kwargs)
            axes_array[-1].append(new_axes)
            if i_horiz == 0:
                last_y_axes = new_axes
            if i_vert == 0:
                last_x_axes = new_axes
    if colorbar:
        colorbar_axes_array = []
        for i_vert in range(n_vert):
            colorbar_axes_array.append([])
            for i_horiz in range(n_horiz):
                x_left = left_margin + (i_horiz + 1) * axis_width \
                    + i_horiz * horiz_sep + colorbar_offset
                y_bottom = bottom_margin \
                    + i_vert * (axis_height + vert_sep) \
                    + axis_height * (1.0 - colorbar_height_fraction) / 2.0
                new_axes = figure.add_axes([x_left / figure_width,
                                            y_bottom / figure_height,
                                            colorbar_width / figure_width,
                                            axis_height * colorbar_height_fraction / figure_height])
                colorbar_axes_array[-1].append(new_axes)
        return (figure, axes_array, colorbar_axes_array)
    else:
        return (figure, axes_array)

def remove_fig_array_axes(axes_array):
    for (i_row, row) in enumerate(axes_array):
        for (i_col, axes) in enumerate(row):
            if i_row > 0:
                for t in axes.get_xticklabels():
                    t.set_visible(False)
            if i_col > 0:
                for t in axes.get_yticklabels():
                    t.set_visible(False)

def find_nearest_index(data, value):
    """Find the index of the entry in data that is closest to value.

    Example:
    >>> data = [0, 3, 5, -2]
    >>> i = partmc.find_nearest_index(data, 3.4)
    returns i = 1

    """
    min_diff = abs(value - data[0])
    min_i = 0
    for i in range(1,len(data)):
        diff = abs(value - data[i])
        if diff < min_diff:
            min_diff = diff
            min_i = i
    return min_i

def label_plot_line(axes, x_data, y_data, label_x, label_text,
                        offset_x=None,
                        offset_y=None,
                        default_offset_mag=5,
                        horizontalalignment="left",
                        verticalalignment="bottom",
                        flip_xy = False,
                        draw_text = True, draw_box = True):
    if offset_x is None:
        if horizontalalignment == "left":
            offset_x = default_offset_mag
        if horizontalalignment == "center":
            offset_x = 0
        if horizontalalignment == "right":
            offset_x = -default_offset_mag
    if offset_y is None:
        if verticalalignment == "bottom":
            offset_y = default_offset_mag
        if verticalalignment == "center":
            offset_y = 0
        if verticalalignment == "top":
            offset_y = -default_offset_mag
    i = find_nearest_index(x_data, label_x)
    label_x = x_data[i]
    label_y = y_data[i]
    kwargs = {"horizontalalignment": horizontalalignment,
              "verticalalignment": verticalalignment,
              }
    if draw_box:
        kwargs["bbox"] = dict(facecolor='white',
                              edgecolor='white')
    axes.annotate(label_text,
                  xy=(label_x, label_y),
                  xycoords='data',
                  xytext=(offset_x, offset_y),
                  textcoords='offset points',
                  **kwargs)

def axes_broken_y(fig, axes, upper_frac=0.5, break_frac=0.05, ybounds=None,
                  xlabel=None, ylabel=None, enforce_upper_frac=False):
    """
    Replace the current axes with a set of upper and lower axes.

    The new axes will be transparent, with a breakmark drawn between them. They
    share the x-axis.  Returns (upper_axes, lower_axes).

    If ybounds=[ymin_lower, ymax_lower, ymin_upper, ymax_upper] is defined,
    upper_frac will be ignored, and the y-axis bounds will be fixed with the
    specified values.
    """
    def breakmarks(fig, axes, y_min, y_max, xwidth=0.03):
        x1, y1, x2, y2 = axes.get_position().get_points().flatten().tolist()
        segment_height = (y_max - y_min) / 3.
        xoffsets = [0, +xwidth, -xwidth, 0]
        yvalues  = [y_min + (i * segment_height) for i in range(4)]
        # Get color of y-axis
        for loc, spine in axes.spines.iteritems():
            if loc  == 'left':
                color = spine.get_edgecolor()
        for x_position in [x1, x2]:
            line = matplotlib.lines.Line2D(
                [x_position + offset for offset in xoffsets], yvalues,
                transform=fig.transFigure, clip_on=False,
                color=color,
                linewidth=0.5)
            axes.add_line(line)

    # Readjust upper_frac if ybounds are defined
    if ybounds:
        if len(ybounds) != 4:
            print "len(ybounds) != 4; aborting..."
            return
        ymin1, ymax1, ymin2, ymax2 = [float(value) for value in ybounds]
        if not enforce_upper_frac:
            data_height1, data_height2 = (ymax1 - ymin1), (ymax2 - ymin2)
            upper_frac = data_height2 / (data_height1 + data_height2)
    x1, y1, x2, y2 = axes.get_position().get_points().flatten().tolist()
    width = x2 - x1
    lower_height = (y2 - y1) * ((1 - upper_frac) - 0.5 * break_frac)
    upper_height = (y2 - y1) * (upper_frac - 0.5 * break_frac)
    upper_bottom = (y2 - y1) - upper_height + y1
    lower_axes = fig.add_axes([x1, y1, width, lower_height],
                              axisbg='None')
    upper_axes = fig.add_axes([x1, upper_bottom, width, upper_height],
                              axisbg='None', sharex=lower_axes)
    # Erase the edges between the axes
    for loc, spine in upper_axes.spines.iteritems():
        if loc == 'bottom':
            spine.set_color('none')
    for loc, spine in lower_axes.spines.iteritems():
        if loc == 'top':
            spine.set_color('none')
    upper_axes.get_xaxis().set_ticks_position('top')
    lower_axes.get_xaxis().set_ticks_position('bottom')
    #plt.setp(upper_axes.get_xticklabels(), visible=False)
    for t in upper_axes.get_xticklabels():
        t.set_visible(False)
    breakmarks(fig, upper_axes, y1 + lower_height, upper_bottom)
    # Set ylims if ybounds are defined
    if ybounds:
        lower_axes.set_ylim(ymin1, ymax1)
        upper_axes.set_ylim(ymin2, ymax2)
        lower_axes.set_autoscaley_on(False)
        upper_axes.set_autoscaley_on(False)
        upper_axes.yaxis.get_label().set_position((0, 1 - (0.5 / (upper_frac/(1+break_frac)))))
        lower_axes.yaxis.get_label().set_position((0, 0.5 / ((1 - upper_frac)/(1+break_frac))))
    # Make original axes invisible
    axes.set_xticks([])
    axes.set_yticks([])
    for loc, spine in axes.spines.iteritems():
        spine.set_color('none')
    return upper_axes, lower_axes
