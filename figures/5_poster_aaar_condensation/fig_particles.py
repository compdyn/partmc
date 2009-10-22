#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math, re
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../../tool")
from pmc_data_nc import *
from pmc_pyx import *
import numpy
from config import *

out_prefix = "figs/particles"
data_prefix = "data/particles"

for [i_run, netcdf_pattern] in netcdf_indexed_patterns:
    out_filename = "%s_%d.pdf" % (out_prefix, i_run)
    print out_filename
    data_filename = "%s_%d.txt" % (data_prefix, i_run)

    filename_list = get_filename_list(netcdf_dir, netcdf_pattern)
    filename = filename_list[0]
    ncf = NetCDFFile(filename)
    env_state = env_state_t(ncf)
    ncf.close()

    data = loadtxt(data_filename)
    data[:,0] /= 60.0 # seconds to minutes
    data[:,1:] *= 1e6 # m to um
    
    min_time_min = data[:,0].min()
    max_time_min = data[:,0].max()
    g = graph.graphxy(
        width = graph_width,
        x = graph.axis.linear(min = 0,
                              max = max_time_min - min_time_min,
                              title = r'time (min)',
                              painter = grid_painter),
        y = graph.axis.log(min = 1e-2,
                           max = 1e2,
                           title = r"diameter ($\rm \mu m$)",
                           painter = major_grid_painter))

    for i in range(1,size(data,1)):
        plot_data = zip(data[:,0] - min_time_min, data[:,i])
        param = float(i - 1) / float(size(data,1) - 2)
        attrs = [rainbow_palette.getcolor(param),
                 style.linewidth.normal]
        g.plot(graph.data.points(plot_data, x = 1, y = 2),
               styles = [graph.style.line(lineattrs = attrs)])

        #write_text_inside(g, show_particles[i]["box label"])

    write_time(g, env_state)

    g.writePDFfile(out_filename)
    print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
    print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
