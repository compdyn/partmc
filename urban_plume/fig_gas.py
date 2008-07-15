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

gas_species = [
    {"species": "O3", "plot": "g1", , "label_time": 8, "label_pos": [1, 1]},
    {"species": "NO2", "plot": "g1", "label_time": 8, "label_pos": [1, 1]},
    {"species": "HCHO", "plot": "g1", "label_time": 8, "label_pos": [1, 1]},
    {"species": "HNO3", "plot": "g2", "label_time": 8, "label_pos": [1, 1]},
    {"species": "SO2", "plot": "g2", "label_time": 8, "label_pos": [1, 1]},
    {"species": "NH3", "plot": "g2", "label_time": 8, "label_pos": [1, 1]},
    ]

netcdf_dir = "out"
netcdf_pattern = r"urban_plume_state_0001_([0-9]{8})\.nc"

gas_state_history = read_history(gas_state_t, netcdf_dir, netcdf_pattern)
env_state = read_any(env_state_t, netcdf_dir, netcdf_pattern)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max([time for [time, gas_state] in gas_state_history]) / 60

c = canvas.canvas()

g1 = c.insert(graph.graphxy(
    width = 6.7,
    x = graph.axis.linear(min = 0.,
                          max = max_time_min,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time
                                               = start_time_of_day_min),
                          title = "local standard time (hours:minutes)",
			  painter = grid_painter),
    y = graph.axis.linear(min = 0.,
                          max = 125,
                          title = "gas concentration (ppb)",
			  painter = grid_painter)))
g2 = c.insert(graph.graphxy(
    width = 6.7,
    ypos = g1.height + 0.5,
    x = graph.axis.linkedaxis(g1.axes["x"],
                              painter = graph.axis.painter.linked(gridattrs = [style.linestyle.dotted])),
    y = graph.axis.linear(min = 0.,
                          max = 20,
                          title = "gas concentration (ppb)",
			  painter = grid_painter)))

graphs = {"g1": g1, "g2": g2}
line_counts = {"g1": 0, "g2": 0}
for i in range(len(gas_species)):
    plot_data = []
    for [time, gas_state] in gas_state_history:
        conc = gas_state.concentration_by_species(gas_species[i]["species"])
        x = time = 60.0
        y = conc
        if y > 0.0:
            plot_data.append([x, y])
    if plot_data != []:
        graph_name = gas_species[i]["plot"]
        g = graphs[graph_name]
        g.plot(graph.data.points(plot_data, x = 1, y = 2,
                                 title = tex_species(gas_species[i])),
               styles = [graph.style.line(lineattrs = [line_style_list[i],
                                                       style.linewidth.THick])])
    else:
        print "warning: only zeros for species %s" % gas_species[i]

c.writePDFfile(out_filename)
print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
