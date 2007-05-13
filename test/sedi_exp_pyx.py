#!/sw/bin/python

import sys, os
sys.path.insert(0, os.path.expanduser("~/.python"))
import pyx
from pyx.graph import *

g = graphxy(width = 15,
	    x = axis.log(min = 1e-7, max = 1e-2,
			 title = "radius (m)"),
	    y = axis.linear(min = 0, max = 1.2e9,
			 title = "number density (\#/m$^3$)"),
	    key = key.key(pos = "tr", dist = 0.1))
g.plot([data.file("out/sedi_exp_mc_out_num_avg.d", x = 2, y = 3,
		  title = "Monte Carlo (0 mins)"),
	data.file("out/sedi_exp_mc_out_num_avg.d", x = 2, y = 8,
		  title = "Monte Carlo (5 mins)"),
	data.file("out/sedi_exp_mc_out_num_avg.d", x = 2, y = 13,
		  title = "Monte Carlo (10 mins)")],
       styles = [style.symbol(size=0.1)])
g.plot([data.file("out/sedi_exp_sect_out_num_avg.d", x = 2, y = 3,
		  title = "Sectional (0 mins)"),
	data.file("out/sedi_exp_sect_out_num_avg.d", x = 2, y = 8,
		  title = "Sectional (5 mins)"),
	data.file("out/sedi_exp_sect_out_num_avg.d", x = 2, y = 13,
		  title = "Sectional (10 mins)")],
       styles = [style.line()])
g.writeEPSfile("out/sedi_exp_pyx_num")
g.writePDFfile("out/sedi_exp_pyx_num")

g = graphxy(width = 5,
	    x = axis.log(min = 1e-7, max = 1e-3,
			 title = "radius (m)",
			 painter = axis.painter.regular(gridattrs=[pyx.attr.changelist([pyx.style.linestyle.dotted, None])])),
	    y = axis.linear(min = 0, max = 8e-6,
			    title = "volume density (m$^3$/m$^3$)",
			    parter = axis.parter.linear([2e-6, 1e-6]),
			    painter = axis.painter.regular(gridattrs=[pyx.style.linestyle.dotted])),
	    key = key.key(pos = "tl", hdist = 0.2, dist = 0.1))
g.plot([data.file("out/sedi_exp_sect_out_vol_avg.d", x = 2, y = 3,
		  title = "sect"),
	data.file("out/sedi_exp_sect_out_vol_avg.d", x = 2, y = 12,
		  title = None)],
       styles = [style.line([pyx.style.linestyle.solid])])
g.plot([data.file("out/sedi_exp_mc_out_vol_avg.d", x = 2, y = 3,
		  title = "MC"),
	data.file("out/sedi_exp_mc_out_vol_avg.d", x = 2, y = 12,
		  title = None)],
       styles = [style.symbol(symbol = style.symbol.square, size=0.05)])
(x,y) = g.pos(2e-5, 6e-6)
g.text(x, y, "0 mins")
(x,y) = g.pos(1e-4, 2.4e-6)
g.text(x, y, "9 mins")
g.writeEPSfile("out/sedi_exp_pyx_vol")
g.writePDFfile("out/sedi_exp_pyx_vol")
