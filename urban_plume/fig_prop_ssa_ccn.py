#!/usr/bin/env python
# Copyright (C) 2007, 2008, 2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
sys.path.append(".")
from fig_helper import *

out_prefix = "figs_prop/ssa_ccn"

ssa_array = loadtxt("SSA output.txt", float)
plot_data_ssa_resolved = []
plot_data_ssa_binned = []
for i in range(ssa_array.shape[0]):
    plot_data_ssa_resolved.append([ssa_array[i,0] * 60, ssa_array[i,1]])
    plot_data_ssa_binned.append([ssa_array[i,0] * 60, ssa_array[i,2]])

ccn_array = loadtxt("figs_prop/ccn.txt", float)
plot_data_ccn_resolved = []
plot_data_ccn_binned = []
for i in range(ccn_array.shape[0]):
    if i % 15 == 1:
        plot_data_ccn_resolved.append([ccn_array[i,0], ccn_array[i,1]])
        plot_data_ccn_binned.append([ccn_array[i,0], ccn_array[i,2]])

max_time_min = 24 * 60
start_time_of_day_min = 6 * 60

for use_color in [True, False]:
    c = canvas.canvas()

    g1 = c.insert(graph.graphxy(
        width = 6.8,
        x = graph.axis.linear(min = 0.,
                              max = max_time_min,
                              parter = graph.axis.parter.linear(tickdists
                                                                = [6 * 60, 3 * 60]),
                              texter = time_of_day(base_time
                                                   = start_time_of_day_min,
                                                   last_suffix = r"\ \ \ \ \ \ "),
                              title = "local standard time (LST) (hours:minutes)",
                              painter = grid_painter),
        y = graph.axis.linear(min = 0.6,
                              max = 1.0,
                              title = r"single scattering albedo",
                              painter = grid_painter)))
    g2 = c.insert(graph.graphxy(
        width = 6.8,
        xpos = g1.xpos + g1.width + 0.6,
        x = graph.axis.linear(min = 0.,
                              max = max_time_min,
                              parter = graph.axis.parter.linear(tickdists
                                                                = [6 * 60, 3 * 60]),
                              texter = time_of_day(base_time
                                                   = start_time_of_day_min,
                                                   first_prefix = r"\ \ \ \ \ \ "),
                              title = "local standard time (LST) (hours:minutes)",
                              painter = grid_painter),
        y2 = graph.axis.linear(min = 0.0,
                              max = 100.0,
                              title = r"unscavenged fraction (\%)",
                              painter = grid_painter)))

    g1.doaxes()
    g2.doaxes()

    if use_color:
        style_attr_resolved = color_list[1]
        style_attr_binned = color_list[2]
    else:
        style_attr_resolved = line_style_list[0]
        style_attr_binned = line_style_list[1]

    g1.plot(
        graph.data.points(plot_data_ssa_resolved, x = 1, y = 2),
        styles = [graph.style.line(lineattrs
                                   = [style_attr_resolved,
                                      style.linewidth.Thick])])
    g1.plot(
        graph.data.points(plot_data_ssa_binned, x = 1, y2 = 2),
        styles = [graph.style.line(lineattrs
                                   = [style_attr_binned,
                                      style.linewidth.Thick])])

    g2.plot(
        graph.data.points(plot_data_ccn_resolved, x = 1, y = 2),
        styles = [graph.style.line(lineattrs
                                   = [style_attr_resolved,
                                      style.linewidth.Thick])])
    g2.plot(
        graph.data.points(plot_data_ccn_binned, x = 1, y = 2),
        styles = [graph.style.line(lineattrs
                                   = [style_attr_binned,
                                      style.linewidth.Thick])])

    label_plot_line_boxed(g1, plot_data_ssa_resolved, 11 * 60,
                          "particle-resolved", [1, 1])
    label_plot_line_boxed(g1, plot_data_ssa_binned, 9.5 * 60,
                          "size-resolved", [0, 0])

    label_plot_line_boxed(g2, plot_data_ccn_resolved, 5 * 60,
                          "particle-resolved", [0, 1])
    label_plot_line_boxed(g2, plot_data_ccn_binned, 5 * 60,
                          r"\parbox{1cm}{size-\\ resolved}", [1, 0])

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    c.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
