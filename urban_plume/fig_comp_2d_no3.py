#!/usr/bin/env python
# Copyright (C) 2007-2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *

times_hour = [1, 6, 12, 24]

netcdf_var = "comp_no3"
netcdf_dim = "composition_no3"
y_axis_label = r"$f_{{\rm BC},{\rm NO_3}}$"
filename = "figs/comp_2d_no3.pdf"

min_val = 0.0
max_val = 2.0

v_space = 0.5
h_space = 0.5

graph_width = 6

data = pmc_var(NetCDFFile("out/urban_plume_with_coag_0001.nc"),
	       netcdf_var,
	       [])
#data.write_summary(sys.stdout)

data.reduce([select("unit", "num_den"),
		 sum("aero_species")])
data.scale_dim(netcdf_dim, 100)
data.scale_dim("dry_radius", 2e6)
data.scale_dim("time", 1.0/3600)

c = canvas.canvas()

g21 = c.insert(graph.graphxy(
    width = graph_width,
    x = graph.axis.log(min = 2.e-3,
                       max = 1.e+0,
                       title = r'dry diameter ($\mu$m)',
                       painter = grid_painter),
    y = graph.axis.linear(min = 0,
                          max = 100,
                          title = y_axis_label,
                          texter = graph.axis.texter.decimal(suffix
                                                             = r"\%"),
                          painter = grid_painter)))
g11 = c.insert(graph.graphxy(
    width = graph_width,
    ypos = g21.height + v_space,
    x = graph.axis.linkedaxis(g21.axes["x"],
                              painter = graph.axis.painter.linked(gridattrs = [style.linestyle.dotted])),
    y = graph.axis.linear(min = 0,
                          max = 100,
                          title = y_axis_label,
                          texter = graph.axis.texter.decimal(suffix
                                                             = r"\%"),
                          painter = grid_painter)))
g22 = c.insert(graph.graphxy(
    width = graph_width,
    xpos = g21.width + h_space,
    x = graph.axis.log(min = 2.e-3,
                       max = 1.e+0,
                       title = r'dry diameter ($\mu$m)',
                       painter = grid_painter),
    y = graph.axis.linkedaxis(g21.axes["y"],
                              painter = graph.axis.painter.linked(gridattrs = [style.linestyle.dotted]))))
g12 = c.insert(graph.graphxy(
    width = graph_width,
    xpos = g11.width + h_space,
    ypos = g22.height + v_space,
    x = graph.axis.linkedaxis(g22.axes["x"],
                              painter = graph.axis.painter.linked(gridattrs = [style.linestyle.dotted])),
    y = graph.axis.linkedaxis(g11.axes["y"],
                              painter = graph.axis.painter.linked(gridattrs = [style.linestyle.dotted]))))

def get_plot_data(time_hour):
    data_slice = module_copy.deepcopy(data)
    data_slice.reduce([select("time", time_hour)])
    data_num = module_copy.deepcopy(data_slice)
    data_num.reduce([sum("dry_radius"), sum(netcdf_dim)])
    data_slice.data = data_slice.data / data_num.data
    plot_data = data_slice.data_2d_list(strip_zero = True,
					min = min_val,
					max = max_val)
    return plot_data

g11.plot(graph.data.points(get_plot_data(times_hour[0]),
                         xmin = 1, xmax = 2, ymin = 3, ymax = 4, color = 5),
         styles = [graph.style.rect(gray_palette)])
g12.plot(graph.data.points(get_plot_data(times_hour[1]),
                         xmin = 1, xmax = 2, ymin = 3, ymax = 4, color = 5),
         styles = [graph.style.rect(gray_palette)])
g21.plot(graph.data.points(get_plot_data(times_hour[2]),
                         xmin = 1, xmax = 2, ymin = 3, ymax = 4, color = 5),
         styles = [graph.style.rect(gray_palette)])
g22.plot(graph.data.points(get_plot_data(times_hour[3]),
                         xmin = 1, xmax = 2, ymin = 3, ymax = 4, color = 5),
         styles = [graph.style.rect(gray_palette)])

g11.doaxes()
g12.doaxes()
g21.doaxes()
g22.doaxes()

g11.dodata()
g12.dodata()
g21.dodata()
g22.dodata()

x_vpos = 0.04
y_vpos = 0.88
boxed_text(g11, x_vpos, y_vpos, "%d hour" % times_hour[0])
boxed_text(g12, x_vpos, y_vpos, "%d hours" % times_hour[1])
boxed_text(g21, x_vpos, y_vpos, "%d hours" % times_hour[2])
boxed_text(g22, x_vpos, y_vpos, "%d hours" % times_hour[3])

add_canvas_color_bar(c,
                     min = min_val,
                     max = max_val,
                     title = r"normalized number density",
                     palette = gray_palette)

c.writePDFfile(filename)
print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
