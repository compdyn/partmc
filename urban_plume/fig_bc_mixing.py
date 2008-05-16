#!/usr/bin/env python
# Copyright (C) 2007-2008 Matthew West, Nicole Riemer
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, glob
import copy as module_copy
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
sys.path.append(os.path.expanduser("~/.python"))
from pyx import *
from Scientific.IO.NetCDF import *

time_hour = 24
composition_lower = [0,  2, 80]
composition_upper = [2, 10, 90]

data_no = pmc_var(NetCDFFile("out/urban_plume_no_coag_0001.nc"),
	       "comp_bc",
	       [])
data_wc = pmc_var(NetCDFFile("out/urban_plume_with_coag_0001.nc"),
	       "comp_bc",
	       [])

#data_no.write_summary(sys.stdout)
#data_wc.write_summary(sys.stdout)

data_no.scale_dim("composition_bc", 100)
data_wc.scale_dim("composition_bc", 100)

data_no.scale_dim("dry_radius", 2e6)
data_wc.scale_dim("dry_radius", 2e6)

data_no.scale_dim("time", 1.0/3600)
data_wc.scale_dim("time", 1.0/3600)

data_no.reduce([select("unit", "mass_den"),
                select("time", time_hour),
		sum("aero_species", without = ["H2O"])])
data_wc.reduce([select("unit", "mass_den"),
                select("time", time_hour),
		sum("aero_species", without = ["H2O"])])

#total_mass = module_copy.deepcopy(data_no)
#total_mass.reduce([sum("dry_radius"), sum("composition_bc")])
#total_mass.scale(1e9)
#print total_mass.data

data_no.scale(1e9) # kg/m^3 to ug/m^3
data_wc.scale(1e9) # kg/m^3 to ug/m^3
data_no.scale(math.log(10.0)) # d/dln(r) to d/dlog10(r)
data_wc.scale(math.log(10.0)) # d/dln(r) to d/dlog10(r)

g = graph.graphxy(
	width = 6.7,
	x = graph.axis.log(title = r'dry diameter ($\mu$m)',
                           min = 0.01, max = 2,
			   painter = grid_painter),
	y = graph.axis.log(title = r'mass density ($\rm \mu g\, m^{-3}$)',
                           min = 1e-6, max = 2e2,
			   painter = grid_painter))
#        key = graph.key.key(pos = "br", vdist = 0.2 * unit.v_cm,
#                            columns = 2))

for i in range(len(composition_lower)):
    data_no_slice = module_copy.deepcopy(data_no)
    data_wc_slice = module_copy.deepcopy(data_wc)
        
    reducer1 = sum("composition_bc", above = composition_lower[i], below = composition_upper[i])
    reducers = [reducer1]
    data_no_slice.reduce(reducers)
    data_wc_slice.reduce(reducers)

    g.plot(graph.data.points(data_no_slice.data_center_list(strip_zero = True),
			   x = 1, y = 2,
                           title = "%d--%d\\%%" % (composition_lower[i], composition_upper[i])),
	   styles = [graph.style.line(lineattrs = [line_style_list[i],style.linewidth.thick])])

    g.plot(graph.data.points(data_wc_slice.data_center_list(strip_zero = True),
			   x = 1, y = 2, 
                           title = "%d--%d\\%%" % (composition_lower[i], composition_upper[i])),
           styles = [graph.style.line(lineattrs = [line_style_list[i],style.linewidth.THick])])

g.doaxes()
g.dodata()

######################################################################
# key

dx = 1 * unit.v_cm
dy = 0.5 * unit.v_cm
xoff = 0.7 * unit.v_cm
yoff = 0.3 * unit.v_cm
length = 0.6 * unit.v_cm
extra = 0.2 * unit.v_cm
coag_off = 0.4 * unit.v_cm
cornerx = 0.3 * unit.v_cm
cornery = 0.3 * unit.v_cm

c = canvas.canvas()

for i in range(len(composition_lower)):
    p = c.text(0, - i * dy,
               "%d--%d\\%% soot" % (composition_lower[i], composition_upper[i]),
               [text.halign.right, text.valign.middle])
    left = p.bbox().left()
    bottom = p.bbox().bottom()
    if i == 0:
        first_top = p.bbox().top()

with_text = c.text(xoff, yoff,
                   "with", [text.halign.center])
without_text = c.text(xoff + dx, yoff,
                      "without", [text.halign.center])
centerx = (with_text.bbox().left() + without_text.bbox().right()) / 2.0
coag_text = c.text(centerx, yoff + coag_off,
                   "coagulation", [text.halign.center])

for i in range(len(composition_lower)):
    c.stroke(path.line(xoff - length / 2.0, - i * dy,
                       xoff + length / 2.0, - i * dy),
                       [line_style_list[i], style.linewidth.THick])
    c.stroke(path.line(xoff + dx - length / 2.0, - i * dy,
                       xoff + dx + length / 2.0, - i * dy),
                       [line_style_list[i], style.linewidth.thick])

#box = path.rect(left - extra, bottom - extra,
#                without_text.bbox().right() - left + 2.0 * extra,
#                coag_text.bbox().top() - bottom + 2.0 * extra)
box = path.path(path.moveto(left - extra, bottom - extra),
                path.lineto(left - extra, first_top + extra),
                path.lineto(with_text.bbox().left() - extra,
                            first_top + extra),
                path.lineto(with_text.bbox().left() - extra,
                            coag_text.bbox().top() + extra),
                path.lineto(without_text.bbox().right() + extra,
                            coag_text.bbox().top() + extra),
                path.lineto(without_text.bbox().right() + extra,
                            bottom - extra),
                path.closepath())
desiredx = g.width - cornerx
desiredy = cornery
transx = - box.bbox().right() + desiredx
transy = - box.bbox().bottom() + desiredy
g.draw(box, [deco.stroked, deco.filled([color.gray.white]),
             trafo.translate(transx, transy)])
g.insert(c, [trafo.translate(transx, transy)])

######################################################################

g.writePDFfile("figs/bc_mixing.pdf")
print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
