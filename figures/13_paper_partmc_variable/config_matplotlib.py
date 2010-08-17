#!/usr/bin/env python
# Copyright (C) 2007-2010 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import sys, os
import math
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt

import config
sys.path.append(config.partmc_tool_dir)
import partmc

matplotlib.rc('text', usetex = True)
matplotlib.rc('xtick.major', pad = 8)
matplotlib.rc('ytick.major', pad = 8)
matplotlib.rc('xtick', labelsize = 10)
matplotlib.rc('legend', fontsize = 10, borderpad = 0.7, borderaxespad = 1)
matplotlib.rc('font', size = 11, family = "serif",
              serif = ["Computer Modern Roman"])
matplotlib.rc('lines', linewidth = 1.5)
matplotlib.rc('patch', linewidth = 0.5)
matplotlib.rc('axes', linewidth = 0.5)

def make_fig(figure_width = 4.1,
             figure_height = None,
             axis_ratio = (1 + math.sqrt(5)) / 2, # golden ratio
             left_margin = 0.8,
             right_margin = 0.2,
             bottom_margin = 0.5,
             top_margin = 0.2,
             colorbar = False,
             colorbar_width = 0.15,
             colorbar_height_fraction = 0.8,
             colorbar_offset = 0.2,
             kind=None,
             ):
    if kind == "1d":
        right_margin = 0.7
    if kind == "2d":
        colorbar = True
        right_margin = 1
        figure_width = 4.4
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
    i = partmc.find_nearest_index(x_data, label_x)
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

def auto_label(axes, lines, labels):
    if len(lines) != len(labels):
        raise Exception("must have same number of lines and labels")

    offset_dist = 0.2

    allowed_anchors = [("left", "bottom", 1, 1),
                       ("center", "bottom", 0, 1),
                       ("right", "bottom", -1, 1),
                       ("left", "center", 1, 0),
                       ("right", "center", -1, 0),
                       ("left", "top", 1, -1),
                       ("center", "top", 0, -1),
                       ("right", "top", -1, -1),
                       ]
    def calc_text_box(label, anchor):
        pass

    def calc_line_point_dist(line_box, x, y, point_x, point_y):
        pass

    def calc_line_box_dist(line_box, x, y, line_x, line_y):
        pass

    def calc_line_dists(line_box, x, y, lines):
        pass

    for i_line in range(len(lines)):
        (line_x, line_y) = lines[i_line]
        for anchor in allowed_anchors:
            for i_attach in len(line_x):
                line_dists = calc_line_dists(line_box, line_x[i_attach],
                                             line_y[i_attach], lines)
                self_dist = line_dists[i_line]
                min_other_dist = min([line_dists[i]
                                      for i in range(len(lines))
                                      if i != i_line])
                
def set_plumetime_axis(axes):
    axes.set_xlim(0, 48)
    axes.set_xticks([0, 12, 24, 36, 48])
    axes.set_xticklabels(["06:00", "18:00", "06:00", "18:00", "06:00"])
    axes.set_xlabel(r"plume time of day (hours:minutes)")

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
    #print upper_axes.yaxis.get_label().get_position()
    #print lower_axes.yaxis.get_label().get_position()
    #print axes.yaxis.get_label().get_position()
    #print axes.yaxis.labelpad
    for loc, spine in axes.spines.iteritems():
        spine.set_color('none')
    return upper_axes, lower_axes

def prepare_efficiency(axes, lower_bound=0.69):
    """
    Set up an efficiency figure with breakmarks to indicate a suppressed zero.

    The y-axis limits are set to (lower_bound, 1.0), as appropriate for an
    efficiency plot, and autoscaling is turned off.
    """
    upper_axes, lower_axes = axes_broken_y(axes, upper_frac=0.97)
    lower_axes.set_yticks([])
    upper_axes.set_ylim(lower_bound, 1.)
    upper_axes.set_autoscaley_on(False)
    return upper_axes, lower_axes

# test these
#ax = plt.axes()
#upper, lower = axes_broken_y(ax, ybounds=[-2., 2.9, 22.1, 30.])
#upper.plot(range(30), range(30))
#lower.plot(range(30), range(30))
#upper.set_ylabel('Data')
#plt.savefig('test')
