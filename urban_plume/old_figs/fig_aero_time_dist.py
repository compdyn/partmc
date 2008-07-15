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

aero_species_1 = [["NO3"], ["NH4"], ["OC"]]
aero_species_2 = [["SO4"], ["BC"], ["ARO1", "ARO2", "ALK1", "OLE1"]]


data = pmc_var(NetCDFFile("out/urban_plume_with_coag_0001.nc"),
	       "aero",
	       [sum("dry_radius"),
		select("unit", "mass_den")])

data.scale_dim("time", 1.0/60)
data.scale(1e9)

c = canvas.canvas()

g2 = c.insert(graph.graphxy(
    width = 6.7,
    x = graph.axis.linear(min = 0.,
                          max = 1440,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time = 6 * 60),
                          title = "local standard time (hours:minutes)",
			  painter = grid_painter),
    y = graph.axis.linear(min = 0.,
                          max = 10,
                          title = r"mass density ($\rm \mu g \, m^{-3}$)",
			  painter = grid_painter)))
g1 = c.insert(graph.graphxy(
    width = 6.7,
    ypos = g2.height + 0.5,
    x = graph.axis.linkedaxis(g2.axes["x"],
                              painter = graph.axis.painter.linked(gridattrs = [style.linestyle.dotted])),
    y = graph.axis.linear(min = 0.,
                          max = 60,
                          title = r"mass density ($\rm \mu g \, m^{-3}$)",
			  painter = grid_painter)))

for i in range(len(aero_species_1)):
    data_slice = module_copy.deepcopy(data)
    data_slice.reduce([sum("aero_species", only = aero_species_1[i])])
    g1.plot(graph.data.points(data_slice.data_center_list(),
			   x = 1, y = 2,
			   title = tex_species(aero_species_1[i])),
	   styles = [graph.style.line(lineattrs = [line_style_list[i], style.linewidth.Thick])])

for i in range(len(aero_species_2)):
    data_slice = module_copy.deepcopy(data)
    data_slice.reduce([sum("aero_species", only = aero_species_2[i])])
    g2.plot(graph.data.points(data_slice.data_center_list(),
			   x = 1, y = 2,
			   title = tex_species(aero_species_2[i])),
	   styles = [graph.style.line(lineattrs = [line_style_list[i], style.linewidth.Thick])])

g1.text(g1.xpos + 1.8,
        g1.ypos + 3.5,
        r"$\rm NO_3$",
        [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g1.text(g1.xpos + 1.9,
        g1.ypos + 1.3,
        r"$\rm NH_4$",
        [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g1.text(g1.xpos + 2.6,
        g1.ypos + 0.7,
        r"POM",
        [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])

g2.text(g2.xpos + 1.6,
        g2.ypos + 3.1,
        r"$\rm SO_4$",
        [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g2.text(g2.xpos + 1.9,
        g2.ypos + 2,
        r"$\rm BC$",
        [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])
g2.text(g2.xpos + 2.4,
        g2.ypos + 1.3,
        r"$\rm SOA$",
        [text.halign.boxleft,text.valign.bottom,color.rgb(0,0,0)])

c.writePDFfile("figs/aero_time_dist.pdf")
print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
