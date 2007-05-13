#!/sw/bin/python

import sys, os
sys.path.insert(0, os.path.expanduser("~/.python"))
from pyx.graph import *

g = graphxy(width = 15,
	    x = axis.log(min = 1e-8, max = 1e-3,
			 title = "radius (m)"),
	    y = axis.log(min = 1e3, max = 1e10,
			 title = "number density (\#/m$^3$)"),
	    key = key.key(pos = "tr", dist = 0.1))
g.plot([data.file("dust_salt_part1_out_num_avg.d", x = 2, y = 3,
		  title = "Monte Carlo (0 s)"),
	data.file("dust_salt_part1_out_num_avg.d", x = 2, y = 7,
		  title = "Monte Carlo (400 s)"),
	data.file("dust_salt_part1_out_num_avg.d", x = 2, y = 11,
		  title = "Monte Carlo (800 s)"),
	data.file("dust_salt_part1_out_num_avg.d", x = 2, y = 13,
		  title = "Monte Carlo (1000 s)")],
       styles = [style.line()])
g.writeEPSfile("dust_salt_num")
g.writePDFfile("dust_salt_num")

g = graphxy(width = 15,
	    x = axis.log(min = 1e-8, max = 1e-3,
			 title = "radius (m)"),
	    y = axis.log(min = 1e-14, max = 1e-4,
			 title = "volume density (m$^3$/m$^3$)"),
	    key = key.key(pos = "tl", dist = 0.1))
g.plot([data.file("dust_salt_part1_out_vol_avg.d", x = 2, y = 3,
		  title = "Monte Carlo (0 s)"),
	data.file("dust_salt_part1_out_vol_avg.d", x = 2, y = 7,
		  title = "Monte Carlo (400 s)"),
	data.file("dust_salt_part1_out_vol_avg.d", x = 2, y = 11,
		  title = "Monte Carlo (800 s)"),
	data.file("dust_salt_part1_out_vol_avg.d", x = 2, y = 13,
		  title = "Monte Carlo (1000s)")],
       styles = [style.line()])
g.writeEPSfile("dust_salt_vol")
g.writePDFfile("dust_salt_vol")

g = graphxy(width = 15,
	    x = axis.log(min = 1e-8, max = 1e-3,
			 title = "radius (m)"),
	    y = axis.log(min = 1e3, max = 1e10,
			 title = "number density (\#/m$^3$)"),
	    key = key.key(pos = "tl", dist = 0.1))
g.plot([data.file("state_dust_salt_part1_0001_00000000_moments_comp.d",
		  x = 2, y = 3,
		  title = "species 1 (0 s)"),
	data.file("state_dust_salt_part1_0001_00000000_moments_comp.d",
		  x = 2, y = 4,
		  title = "species 2 (0 s)"),
	data.file("state_dust_salt_part2_0001_00000800_moments_comp.d",
		  x = 2, y = 3,
		  title = "species 1 (800 s)"),
	data.file("state_dust_salt_part2_0001_00000800_moments_comp.d",
		  x = 2, y = 4,
		  title = "species 2 (800 s)"),
	data.file("state_dust_salt_part2_0001_00000800_moments_comp.d",
		  x = 2, y = 5,
		  title = "mixed (800 s)"),
	data.file("state_dust_salt_part2_0001_00001000_moments_comp.d",
		  x = 2, y = 3,
		  title = "species 1 (1000 s)"),
	data.file("state_dust_salt_part2_0001_00001000_moments_comp.d",
		  x = 2, y = 4,
		  title = "species 2 (1000 s)"),
	data.file("state_dust_salt_part2_0001_00001000_moments_comp.d",
		  x = 2, y = 5,
		  title = "mixed (1000s)")],
       styles = [style.line()])
g.writeEPSfile("dust_salt_comp")
g.writePDFfile("dust_salt_comp")
