#!/sw/bin/python

import sys, os
sys.path.insert(0, os.path.expanduser("~/.python"))
import pyx
from pyx.graph import *

g = graphxy(width = 5,
	    x = axis.log(min = 1e-7, max = 1e-3,
			 title = "radius (m)",
			 painter = axis.painter.regular(gridattrs=[pyx.attr.changelist([pyx.style.linestyle.dotted, None])])),
	    y = axis.log(min = 1e4, max = 1e10,
			 title = "number density (\#/m$^3$)",
			 painter = axis.painter.regular(gridattrs=[pyx.style.linestyle.dotted])),
	    key = key.key(pos = "tl", hdist = 0.2, dist = 0.1))
g.plot([data.file("out_golovin_exact_num_avg.d", x = 2, y = 3,
		  title = "analyt"),
	data.file("out_golovin_exact_num_avg.d", x = 2, y = 8, title = None),
	data.file("out_golovin_exact_num_avg.d", x = 2, y = 13, title = None)],
       styles = [style.line(lineattrs=[pyx.style.linestyle.solid])])
g.plot([data.file("out_golovin_mc_num_avg.d", x = 2, y = 3,
		  title = "MC"),
	data.file("out_golovin_mc_num_avg.d", x = 2, y = 8, title = None),
	data.file("out_golovin_mc_num_avg.d", x = 2, y = 13, title = None)],
       styles = [style.symbol(symbol = style.symbol.square, size=0.05)])
(x,y) = g.pos(2e-5, 5e8)
g.text(x, y, "0 mins")
(x,y) = g.pos(4e-5, 5e7)
g.text(x, y, "5 mins")
(x,y) = g.pos(8e-5, 5e6)
g.text(x, y, "10 mins")
g.writeEPSfile("plot_golovin_pyx_num")
g.writePDFfile("plot_golovin_pyx_num")

g = graphxy(width = 6,
	    x = axis.log(min = 1e-7, max = 1e-3,
			 title = "radius (m)",
			 painter = axis.painter.regular(gridattrs=[pyx.style.linestyle.dotted])),
	    y = axis.log(min = 1e-14, max = 1e-4,
			 title = "volume density (m$^3$/m$^3$)",
			 painter = axis.painter.regular(gridattrs=[pyx.style.linestyle.dotted])),
	    key = key.key(pos = "tl", dist = 0.1))
g.plot([data.file("out_golovin_mc_vol_avg.d", x = 2, y = 3,
		  title = "Monte Carlo (0 mins)"),
	data.file("out_golovin_mc_vol_avg.d", x = 2, y = 8,
		  title = "Monte Carlo (5 mins)"),
	data.file("out_golovin_mc_vol_avg.d", x = 2, y = 13,
		  title = "Monte Carlo (10 mins)")],
       styles = [style.symbol(size=0.1)])
g.plot([data.file("out_golovin_exact_vol_avg.d", x = 2, y = 3,
		  title = "Analytical (0 mins)"),
	data.file("out_golovin_exact_vol_avg.d", x = 2, y = 8,
		  title = "Analytical (5 mins)"),
	data.file("out_golovin_exact_vol_avg.d", x = 2, y = 13,
		  title = "Analytical (10 mins)")],
       styles = [style.line()])
g.writeEPSfile("plot_golovin_pyx_vol")
g.writePDFfile("plot_golovin_pyx_vol")
