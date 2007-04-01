#!/sw/bin/python

import sys, os
sys.path.insert(0, os.path.expanduser("~/.python"))
from pyx.graph import *

g = graphxy(width = 15,
	    x = axis.log(min = 1e-7, max = 1e-3,
			 title = "radius (m)"),
	    y = axis.log(min = 1e4, max = 1e10,
			 title = "number density (\#/m$^3$)"),
	    key = key.key(pos = "tr", dist = 0.1))
g.plot([data.file("out_golovin_mc_num_avg.d", x = 2, y = 3,
		  title = "Monte Carlo (0 mins)"),
	data.file("out_golovin_mc_num_avg.d", x = 2, y = 8,
		  title = "Monte Carlo (5 mins)"),
	data.file("out_golovin_mc_num_avg.d", x = 2, y = 13,
		  title = "Monte Carlo (10 mins)")],
       styles = [style.symbol(size=0.1)])
g.plot([data.file("out_golovin_exact_num_avg.d", x = 2, y = 3,
		  title = "Analytical (0 mins)"),
	data.file("out_golovin_exact_num_avg.d", x = 2, y = 8,
		  title = "Analytical (5 mins)"),
	data.file("out_golovin_exact_num_avg.d", x = 2, y = 13,
		  title = "Analytical (10 mins)")],
       styles = [style.line()])
g.writeEPSfile("tpn")
g.writePDFfile("tpn")

g = graphxy(width = 15,
	    x = axis.log(min = 1e-7, max = 1e-3,
			 title = "radius (m)"),
	    y = axis.log(min = 1e-14, max = 1e-4,
			 title = "volume density (m$^3$/m$^3$)"),
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
g.writeEPSfile("tpv")
g.writePDFfile("tpv")
