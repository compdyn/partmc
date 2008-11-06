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
from fig_helper import *

out_prefix = "figs/test_emission"

emission_netcdf_dir = r"../test/emission/out"

mc_data = pmc_var(NetCDFFile(os.path.join(emission_netcdf_dir,
                                          "emission_mc_0001.nc")),
		  "aero",
		  [sum("aero_species")])
exact_data = pmc_var(NetCDFFile(os.path.join(emission_netcdf_dir,
                                             "emission_exact_0001.nc")),
		     "aero",
		     [sum("aero_species")])
sect_data = pmc_var(NetCDFFile(os.path.join(emission_netcdf_dir,
                                            "emission_sect_0001.nc")),
		    "aero",
		    [sum("aero_species")])

######################################################################

samples = [
    {"radius": 1.5e-5,
     "label": "initial", "label_time": 2, "label_pos": [0, 1],
     "label_offset": unit.v_mm},
    {"radius": 3.0e-5,
     "label": "emissions", "label_time": 22, "label_pos": [1, 1],
     "label_offset": 1.5 * unit.v_mm},
    {"radius": 5.0e-5,
     "label": "background", "label_time": 17, "label_pos": [1, 1],
     "label_offset": unit.v_mm},
    ]

for use_color in [True, False]:
    g = graph.graphxy(width = 6,
                      x = graph.axis.linear(min = 0,
                                            max = 24,
                                            parter = graph.axis.parter.linear(tickdists = [6, 3]),
                                            title = r"elapsed time (hours)",
                                            painter = grid_painter),
                      y = graph.axis.linear(min = 0,
                                            max = 1.5e10,
                                            title = r"number concentration ($\rm m^{-3}$)",
                                            painter = grid_painter))

    for i in range(len(samples)):
        data_slice = module_copy.deepcopy(mc_data)
        data_slice.reduce([select("unit", "num_den"),
                           select("radius", samples[i]["radius"])])
        data_slice.scale_dim("time", 1.0/3600)
        plot_data = data_slice.data_center_list()
        plot_data = [plot_data[k] for k in range(len(plot_data))
                     if k % 70 == 0]
        if use_color:
            plot_color = color_list[i]
        else:
            plot_color = color.gray.black
        g.plot(graph.data.points(plot_data, x = 1, y = 2),
               styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
                                            size = 0.08,
                                            symbolattrs = [plot_color])])
        ##################
        data_slice = module_copy.deepcopy(exact_data)
        data_slice.reduce([select("unit", "num_den"),
                           select("radius", samples[i]["radius"])])
        data_slice.scale_dim("time", 1.0/3600)
        plot_data = data_slice.data_center_list()
        g.plot(graph.data.points(plot_data, x = 1, y = 2),
               styles = [graph.style.line(lineattrs = [plot_color,
                                                       style.linestyle.solid])])
        label_plot_line(g, plot_data, samples[i]["label_time"],
                        samples[i]["label"], label_pos = samples[i]["label_pos"],
                        label_offset = samples[i]["label_offset"])

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    g.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
