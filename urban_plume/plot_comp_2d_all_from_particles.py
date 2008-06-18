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

f = NetCDFFile("out/urban_plume_with_coag_state_0001_00000024.nc")
particles = read_particles(f)
x_axis = pmc_log_axis(min = 1e-2, max = 1, n_bin = 160)
y_axis = pmc_linear_axis(min = 0, max = 1, n_bin = 100)
bin_array = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
a_species = ["SO4", "NO3", "Cl", "NH4", "MSA", "ARO1", "ARO2", "ALK1",
             "OLE1", "API1", "API2", "LIM1", "LIM2", "CO3", "Na", "Ca",
             "OIN", "OC"]
b_species = ["BC"]
for particle in particles:
    diameter = 2.0 * particle.radius() * 1e6 # um
    comp_frac = particle.comp_frac(a_species, b_species, "mass") * 100
    x_bin = x_axis.find(diameter)
    y_bin = y_axis.find(comp_frac)
    bin_array[x_bin, y_bin] += 1.0 / particle.comp_vol \
                               / x_axis.grid_size(x_bin) \
                               / y_axis.grid_size(y_bin)
bin_array = bin_array / bin_array.max()

g = graph.graphxy(
    width = 5,
    x = graph.axis.log(min = x_axis.min,
                       max = x_axis.max,
                       title = r'diameter ($\mu$m)'),
    y = graph.axis.linear(min = y_axis.min,
                          max = y_axis.max,
                          title = r"$f_{\rm BC,all}$",
                          texter = graph.axis.texter.decimal(suffix = r"\%")))

g.plot(graph.data.points(pmc_histogram_2d(bin_array, x_axis, y_axis),
                         xmin = 1, xmax = 2, ymin = 3, ymax = 4, color = 5),
       styles = [graph.style.rect(gray_palette)])
g.writePDFfile("out/comp_2d_all_from_particles.pdf")
