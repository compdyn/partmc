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
sys.path.append(".")
from fig_helper import *

gas_species = [
    {"species": "O3", "plot": "g1", "label_time": 6.1, "label_pos": [1, 1], "label_offset": 0 * unit.v_mm},
    {"species": "NO2", "plot": "g1", "label_time": 8, "label_pos": [1, 1], "label_offset": 0.5 * unit.v_mm},
    {"species": "HCHO", "plot": "g2", "label_time": 8, "label_pos": [0, 0], "label_offset": 0 * unit.v_mm},
    {"species": "HNO3", "plot": "g2", "label_time": 8, "label_pos": [1, 1], "label_offset": 0 * unit.v_mm},
    {"species": "SO2", "plot": "g2", "label_time": 15, "label_pos": [0, 1], "label_offset": 0 * unit.v_mm},
    {"species": "NH3", "plot": "g2", "label_time": 11, "label_pos": [0, 1], "label_offset": 0 * unit.v_mm},
    ]

out_prefix = "figs/gas"

gas_state_history = read_history(gas_state_t, netcdf_dir_wc, netcdf_pattern_wc)
env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max([time for [time, gas_state] in gas_state_history]) / 60

for use_color in [True, False]:
    c = canvas.canvas()

    g2 = c.insert(graph.graphxy(
        width = 6.7,
        x = graph.axis.linear(min = 0.,
                              max = max_time_min,
                              parter = graph.axis.parter.linear(tickdists
                                                                = [6 * 60, 3 * 60]),
                              texter = time_of_day(base_time
                                                   = start_time_of_day_min),
                              title = "local standard time (LST) (hours:minutes)",
                              painter = grid_painter),
        y = graph.axis.linear(min = 0.,
                              max = 20,
                              title = "gas mole fraction (ppb)",
                              painter = grid_painter)))
    g1 = c.insert(graph.graphxy(
        width = 6.7,
        ypos = g2.height + 0.5,
        x = graph.axis.linkedaxis(g2.axes["x"],
                                  painter = graph.axis.painter.linked(gridattrs = [style.linestyle.dotted])),
        y = graph.axis.linear(min = 0.,
                              max = 150,
                              title = "gas mole fraction (ppb)",
                              painter = grid_painter)))

    g1.doaxes()
    g2.doaxes()

    graphs = {"g1": g1, "g2": g2}
    line_counts = {"g1": 0, "g2": 0}
    for i in range(len(gas_species)):
        plot_data = []
        for [time, gas_state] in gas_state_history:
            conc = gas_state.concentration_by_species(gas_species[i]["species"])
            if conc > 0.0:
                plot_data.append([time / 60.0, conc])
        if plot_data != []:
            graph_name = gas_species[i]["plot"]
            g = graphs[graph_name]
            if use_color:
                style_attr = color_list[line_counts[graph_name]]
            else:
                style_attr = line_style_list[line_counts[graph_name]]
            g.plot(graph.data.points(plot_data, x = 1, y = 2,
                                     title = tex_species(gas_species[i])),
                   styles = [graph.style.line(
                lineattrs = [style_attr,
                             style.linewidth.THick])])
            line_counts[graph_name] += 1
            label = tex_species(gas_species[i]["species"])
            label_plot_line_boxed(g, plot_data, gas_species[i]["label_time"] * 60.0,
                        label, gas_species[i]["label_pos"], label_offset = gas_species[i]["label_offset"])
        else:
            print "warning: only zeros for species %s" % gas_species[i]

        a_data = array([v for [t,v] in plot_data])
        j = a_data.argmax()
        if not use_color:
            print "%s max = %g ppb at %s LST" \
                  % (gas_species[i]["species"], plot_data[j][1],
                     time_of_day_string(plot_data[j][0] * 60
                                        + env_state.start_time_of_day))

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    c.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
