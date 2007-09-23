#!/usr/bin/env python
# Copyright (C) 2007 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
sys.path.append("../tool")
from pmc_data import *
from pmc_pyx import *
sys.path.append(os.path.expanduser("~/.python"))
from pyx import *

######################################################################

g = graph.graphxy(width = 10,
		  x = graph.axis.linear(title = "time (s)",
					painter = grid_painter),
		  y = graph.axis.linear(title = "number density of small particles (\#/m$^3$)",
					painter = grid_painter),
		  key = graph.key.key(pos = "tr"))
g.plot(graph.data.file("out/sedi_bidisperse_mc_counts.d",
		       x = 1, y = 2,
		       title = "Monte Carlo"),
       styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				    size = 0.05,
				    symbolattrs = [color_list[0]])])
g.plot(graph.data.file("out/sedi_bidisperse_ode_counts.d",
		       x = 1, y = 2,
		       title = "ODE"),
       styles = [graph.style.line(lineattrs = [color_list[1]])])
g.writePDFfile("out/sedi_bidisperse_num_small")

######################################################################

file = open("out/sedi_bidisperse_mc_counts.d")
data = []
for line in file:
    data.append([float(x) for x in line.split()])
data_small = [[d[0], d[1]] for d in data if d[1] > 0]
file.close()

g = graph.graphxy(width = 10,
		  x = graph.axis.linear(title = "time (s)",
					painter = grid_painter),
		  y = graph.axis.log(title = "number density of small particles (\#/m$^3$)",
				     painter = grid_painter),
		  key = graph.key.key(pos = "tr"))
g.plot(graph.data.list(data_small, x = 1, y = 2,
		       title = "Monte Carlo"),
       styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				    size = 0.05,
				    symbolattrs = [color_list[0]])])
g.plot(graph.data.file("out/sedi_bidisperse_ode_counts.d",
		       x = 1, y = 2,
		       title = "ODE"),
       styles = [graph.style.line(lineattrs = [color_list[1]])])
g.writePDFfile("out/sedi_bidisperse_num_small_log")

######################################################################

g = graph.graphxy(width = 10,
		  x = graph.axis.linear(title = "time (s)",
					painter = grid_painter),
		  y = graph.axis.linear(title = "volume density of big particles (m$^3$/m$^3$)",
					painter = grid_painter),
		  key = graph.key.key(pos = "br"))
g.plot(graph.data.file("out/sedi_bidisperse_mc_counts.d",
		       x = 1, y = 3,
		       title = "Monte Carlo"),
       styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
				    size = 0.05,
				    symbolattrs = [color_list[0]])])
g.plot(graph.data.file("out/sedi_bidisperse_ode_counts.d",
		       x = 1, y = 3,
		       title = "ODE"),
       styles = [graph.style.line(lineattrs = [color_list[1]])])
g.writePDFfile("out/sedi_bidisperse_vol_big")
