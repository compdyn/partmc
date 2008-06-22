#!/usr/bin/env python
# Copyright (C) 2007-2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
import numpy

graph_width = 8
v_space = 0.5

f = NetCDFFile("out/runs_20080620/data/20080620_gas_halved_3am/urban_plume_0.5_3am_state_0001_00000024.nc")
particles = read_particles(f)
x_axis = pmc_log_axis(min = 1e-2, max = 2, n_bin = 160)
y_axis = pmc_linear_axis(min = 0, max = 100, n_bin = 100)
dry_bin_array = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
wet_bin_array = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
a_species = ["SO4", "NO3", "Cl", "NH4", "MSA", "ARO1", "ARO2", "ALK1",
             "OLE1", "API1", "API2", "LIM1", "LIM2", "CO3", "Na", "Ca",
             "OIN", "OC"]
b_species = ["BC"]
for particle in particles:
    diameter = 2.0 * particle.radius() * 1e6 # um
    comp_frac = particle.comp_frac(a_species, b_species, "mass") * 100
    x_bin = x_axis.find(diameter)
    y_bin = y_axis.find(comp_frac)
    if particle.species_mass("H2O") == 0.0:
        dry_bin_array[x_bin, y_bin] = 1.0
    else:
        wet_bin_array[x_bin, y_bin] = particle.species_mass("H2O") \
                                      / particle.mass()
max_wet = wet_bin_array.max()
wet_bin_array = wet_bin_array / max_wet

c = canvas.canvas()

g2 = c.insert(graph.graphxy(
    width = graph_width,
    x = graph.axis.log(min = x_axis.min,
                       max = x_axis.max,
                       title = r'diameter ($\mu$m)'),
    y = graph.axis.linear(min = y_axis.min,
                          max = y_axis.max,
                          title = r"$f_{\rm BC,all}$",
                          texter = graph.axis.texter.decimal(suffix = r"\%"))))
g1 = c.insert(graph.graphxy(
    width = graph_width,
    ypos = g2.height + v_space,
    x = graph.axis.linkedaxis(g2.axes["x"]),
    y = graph.axis.linear(min = y_axis.min,
                          max = y_axis.max,
                          title = r"$f_{\rm BC,all}$",
                          texter = graph.axis.texter.decimal(suffix = r"\%"))))

g1.plot(graph.data.points(pmc_histogram_2d(dry_bin_array, x_axis, y_axis),
                         xmin = 1, xmax = 2, ymin = 3, ymax = 4, color = 5),
       styles = [graph.style.rect(rainbow_palette)])
g2.plot(graph.data.points(pmc_histogram_2d(wet_bin_array, x_axis, y_axis),
                         xmin = 1, xmax = 2, ymin = 3, ymax = 4, color = 5),
       styles = [graph.style.rect(rainbow_palette)])
add_color_bar(g2, min = 0.0, max = max_wet * 100,
              title = r"water fraction", palette = rainbow_palette,
              bar_x_offset = 0.8,
              texter = graph.axis.texter.decimal(suffix = r"\%"))

boxed_text(g1, 0.8, 0.9, "dry")
boxed_text(g2, 0.8, 0.9, "wet")

c.writePDFfile("out/comp_2d_all_water.pdf")
