#!/usr/bin/env python
# Copyright (C) 2007, 2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *

env = ["H"]

subdir = "run_RH"
if len(sys.argv) > 1:
    subdir = sys.argv[1]

data = pmc_var(NetCDFFile("out/%s/urban_plume_RH05_0001.nc" % subdir),
	       "env_state",
	       [])

data.write_summary(sys.stdout)

data.scale_dim("time", 1.0/60)

temp_data = module_copy.deepcopy(data)
temp_data.reduce([select("env", "temp")])
rh_data = module_copy.deepcopy(data)
rh_data.reduce([select("env", "rel_humid")])
rh_data.scale(100.0)
print rh_data.data
height_data = module_copy.deepcopy(data)
height_data.reduce([select("env", "height")])

g = graph.graphxy(
    width = 10,
    x = graph.axis.linear(min = 0.,
                          max = 720,
#                          parter = graph.axis.parter.linear(tickdists = [6, 3]),
			  title = "time (hour)",
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time = 6 * 60),
			  painter = grid_painter),
    y = graph.axis.linear(title = "temperature (K)",
                          painter = grid_painter),
    y2 = graph.axis.linear(title = r"relative humidity (\%)"),
    y4 = graph.axis.linear(title = "mixing height (m)"),
    key = graph.key.key(pos = "tr"))

g.plot(graph.data.list(temp_data.data_center_list(),
			   x = 1, y = 2,
                           title = "temperature"),
             styles = [graph.style.line(lineattrs = [color_list[1], style.linewidth.THick])])

g.plot(graph.data.list(rh_data.data_center_list(),
			   x = 1, y2 = 2,
                           title = "relative humidity"),
             styles = [graph.style.line(lineattrs = [color_list[2],style.linewidth.THick])])

g.plot(graph.data.list(height_data.data_center_list(),
			   x = 1, y4 = 2,
                           title = "mixing height"),
             styles = [graph.style.line(lineattrs = [color_list[5],style.linewidth.THick])])

g.writePDFfile("out/%s/figs_RH05/env_RH05.pdf" % subdir)
