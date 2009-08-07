#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import random as py_random
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append(".")
from fig_helper import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *

out_prefix = "figs_aging/aging_aero_2d_aging_time"

max_val = 4.0

const = load_constants("../src/constants.f90")

coag_suffix = "wc"
bin = level_mid + 1

filename = os.path.join(aging_data_dir,
                        "particle_aging_%s_plot_data_%08d.txt" % (coag_suffix, bin))
value = loadtxt(filename)
value = value / max_val
value = value.clip(0.0, 1.0)

plot_data = pmc_histogram_2d(value, diameter_axis, aging_time_axis)

for color in [True, False]:
    g = graph.graphxy(
        width = grid_graph_width,
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
                           title = diameter_axis_label),
        y = graph.axis.log(min = aging_time_axis.min,
                              max = aging_time_axis.max,
                              title = r"aging time $t_{\rm aging}\ (\rm h)$"))

    if color:
        palette = rainbow_palette
    else:
        palette = nonlinear_gray_palette

    g.dolayout()

    for axisname in ["x", "y"]:
        for t in g.axes[axisname].data.ticks:
            if t.ticklevel is not None:
                g.stroke(g.axes[axisname].positioner.vgridpath(t.temp_v),
                         [style.linestyle.dotted])
                    
    g.plot(graph.data.points(plot_data,
                             xmin = 1, xmax = 2, ymin = 3, ymax = 4,
                             color = 5),
           styles = [hsb_rect(palette)])

    g.dodata()
    g.doaxes()

    add_canvas_color_bar(
        g,
        min = 0.0,
        max = max_val,
        xpos = g.xpos + g.width + grid_h_space,
        ybottom = g.ypos,
        ytop = g.ypos + g.height,
        title = r"normalized number conc.",
        palette = palette)

    if color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    g.writePDFfile(out_filename)
    if not color:
        print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
