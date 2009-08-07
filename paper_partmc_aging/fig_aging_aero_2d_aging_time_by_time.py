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

const = load_constants("../src/constants.f90")

coag_suffix = "wc"
bin = level_mid + 1

for time_type in ["emission", "aging"]:
    filename = os.path.join(aging_data_dir,
                            "particle_aging_by_time_%s_%s_plot_data_%08d.txt" % (time_type, coag_suffix, bin))

    data = loadtxt(filename)
    diameter = data[:,0]
    aging_time = data[:,1]
    z_time = data[:,2]

    z_time_max = z_time.max()
    z_time_scaled = z_time / z_time_max

    plot_data = grid_plot_data(diameter, aging_time, z_time_scaled,
                               diameter_axis, aging_time_axis,
                               diameter_axis_hi, aging_time_axis_hi)

    for color in [True, False]:
        g = graph.graphxy(
            width = grid_graph_width,
            x = graph.axis.log(min = diameter_axis_min,
                               max = diameter_axis_max,
                               title = r"dry diameter at emission $D\ (\rm\mu m)$"),
            y = graph.axis.linear(min = aging_time_axis.min,
                                  max = aging_time_axis.max,
                                  parter = graph.axis.parter.linear(tickdists = [6, 3]),
                                  title = r"aging time $t_{\rm aging}\ (\rm h)$"))

        if color:
            palette = rainbow_palette
        else:
            palette = gray_palette

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

        write_text_outside(g, r"critical supersaturation $S_{\rm c} = %.1f\%%$" % (level_mid_value * 100))

        add_canvas_color_bar(
            g,
            min = 0,
            max = z_time_max / 3600.0,
            xpos = g.xpos + g.width + grid_h_space,
            ybottom = g.ypos,
            ytop = g.ypos + g.height,
            title = r"elapsed time at %s $(h)$" % time_type,
            palette = palette)

        if color:
            out_filename = "%s_%s_%s_color.pdf" % (out_prefix, time_type, coag_suffix)
        else:
            out_filename = "%s_%s_%s_bw.pdf" % (out_prefix, time_type, coag_suffix)
        g.writePDFfile(out_filename)
        if not color:
            print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
            print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
