#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../../tool")
from pmc_data_nc import *
from pmc_pyx import *
from config import *

out_filename_max_ss = "figs/time_max_ss.pdf"
out_filename_active = "figs/time_active.pdf"

particle_data = [
    [1,  0.1741616308,  41.6524],
    [7,  0.09059989756, 39.5619],
    [15, 0.07041647553, 43.5748],
    [24, 0.09695690379, 56.1   ],
    ]

size_avg_data = [
    [1,  0.1525054598,  47.6125],
    [7,  0.08602255906, 37.1262],
    [15, 0.04186271645, 58.8006],
    [24, 0.05368977188, 62.2821],
    ]

comp_avg_data = [
    [1,  ],
    [7,  ],
    [15, ],
    [24, ],
    ]
     
g_max_ss = graph.graphxy(
    width = graph_width,
    x = graph.axis.linear(min = 0,
                          max = 24,
                          title = "time (h)",
                          painter = grid_painter,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6, 3])),
    y = graph.axis.linear(min = 0,
                          max = 0.2,
                          title = "maximum supersaturation ($\%$)",
                          painter = grid_painter))

g_active = graph.graphxy(
    width = graph_width,
    x = graph.axis.linear(min = 0,
                          max = 24,
                          title = "time (h)",
                          painter = grid_painter,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6, 3])),
    y = graph.axis.linear(min = 0,
                          max = 100,
                          title = "activated fraction ($\%$)",
                          painter = grid_painter))

g_max_ss.doaxes()
g_active.doaxes()

for (i_data, (data, name)) in enumerate(zip([particle_data, size_avg_data],
                                            ["particle", "size-avg"])):
    max_ss_data = [[d[0], d[1]] for d in data]
    active_data = [[d[0], d[2]] for d in data]
    g_max_ss.plot(graph.data.points(max_ss_data, x = 1, y = 2),
                  styles = [graph.style.line(lineattrs = [color_list[i_data]])])
    g_active.plot(graph.data.points(active_data, x = 1, y = 2),
                  styles = [graph.style.line(lineattrs = [color_list[i_data]])])
    label_plot_line_boxed(g_max_ss, max_ss_data, 15, name, [1, 1])
    label_plot_line_boxed(g_active, active_data, 15, name, [1, 1])

g_max_ss.writePDFfile(out_filename_max_ss)
g_active.writePDFfile(out_filename_active)
