#!/usr/bin/env python
# Copyright (C) 2007 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
sys.path.append("../../tool")
from pmc_pyx import *
sys.path.append(os.path.expanduser("~/.python"))
from pyx import *

n_list = [1, 4, 10, 20, 30]

g = graph.graphxy(width = 10,
		  x = graph.axis.linear(title = "n",
					painter = grid_painter),
		  y = graph.axis.linear(title = "Prob(n)",
					painter = grid_painter),
		  key = graph.key.key(pos = "tr"))

for i in range(len(n_list)):
    g.plot(graph.data.file("out/poisson_%d.d" % n_list[i],
			   x = 1, y = 2,
			   title = "rate %d (exact)" % n_list[i]),
	   styles = [graph.style.line(lineattrs = [color_list[i]])])
    g.plot(graph.data.file("out/poisson_%d.d" % n_list[i],
			   x = 1, y = 3,
			   title = "rate %d (sampled)" % n_list[i]),
	   styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
					size = 0.05,
					symbolattrs = [color_list[i]])])

g.writePDFfile("out/poisson")
