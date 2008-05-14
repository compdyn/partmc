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

gas_species_1 = ["O3", "NO2", "HCHO"]
gas_species_2 = ["HNO3", "SO2", "NH3"]

data = pmc_var(NetCDFFile("out/urban_plume_with_coag_0001.nc"),
	       "gas",
	       [])

data.scale_dim("time", 1.0/60)

c = canvas.canvas()

g1 = c.insert(graph.graphxy(
    width = 6.5,
    x = graph.axis.linear(min = 0.,
                          max = 1440,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time = 6 * 60),
                          title = "local standard time",
			  painter = grid_painter),
    y = graph.axis.linear(min = 0.,
                          max = 125,
                          title = "gas concentration (ppb)",
			  painter = grid_painter)))
g2 = c.insert(graph.graphxy(
    width = 6.5,
    ypos = g1.height + 0.5,
    x = graph.axis.linkedaxis(g1.axes["x"],
                              painter = graph.axis.painter.linked(gridattrs = [style.linestyle.dotted])),
    y = graph.axis.linear(min = 0.,
                          max = 20,
                          title = "gas concentration (ppb)",
			  painter = grid_painter)))

for i in range(len(gas_species_1)):
    data_slice = module_copy.deepcopy(data)
    data_slice.reduce([select("gas_species", gas_species_1[i])])
    g1.plot(graph.data.points(data_slice.data_center_list(),
			   x = 1, y = 2,
			   title = tex_species(gas_species_1[i])),
	   styles = [graph.style.line(lineattrs = [line_style_list[i], style.linewidth.Thick])])

for i in range(len(gas_species_2)):
    data_slice = module_copy.deepcopy(data)
    data_slice.reduce([select("gas_species", gas_species_2[i])])
    g2.plot(graph.data.points(data_slice.data_center_list(),
			   x = 1, y = 2,
			   title = tex_species(gas_species_2[i])),
	   styles = [graph.style.line(lineattrs = [line_style_list[i], style.linewidth.Thick])])

g1.text(g1.xpos + 1.5,
        g1.ypos + 3.5,
        r"$\rm O_3$",
        [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g1.text(g1.xpos + 2,
        g1.ypos + 2.2,
        r"$\rm NO_2$",
        [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g1.text(g1.xpos + 2.6,
        g1.ypos + 0.8,
        r"HCHO",
        [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])

g2.text(g2.xpos + 1.5,
        g2.ypos + 3.1,
        r"$\rm HNO_3$",
        [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g2.text(g2.xpos + 2.2,
        g2.ypos + 2,
        r"$\rm SO_2$",
        [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g2.text(g2.xpos + 2.4,
        g2.ypos + 0.7,
        r"$\rm NH_3$",
        [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])

c.writePDFfile("figs/gas_time_dist.pdf")
print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
