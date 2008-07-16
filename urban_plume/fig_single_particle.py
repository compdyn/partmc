#!/usr/bin/env python
# Copyright (C) 2007-2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math, re
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
import numpy
from fig_helper import *

out_prefix = "figs/particle"

aero_species = [["", ["NO3"]],
                ["", ["NH4"]],
                ["", ["OC"]],
                ["", ["H2O"]],
                ["", ["SO4"]],
                ["", ["BC"]],
                ["SOA", ["ARO1", "ARO2", "ALK1", "OLE1"]],
                ]

particle_ids = [p["id"] for p in show_particles]
particle_history = read_history(lambda ncf:
                                read_particles(ncf, ids = particle_ids),
                                netcdf_dir_wc, netcdf_pattern_wc)
env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max([time for [time, x] in particle_history]) / 60

def particle_by_id(particle_list, id):
    for particle in particle_list:
        if particle.id == id:
            return particle
    return None

for i in range(len(show_particles)):
    g = graph.graphxy(
        width = 6.7,
        x = graph.axis.linear(min = 0,
                              max = max_time_min,
                              title = r'time (LST)',
                              parter = graph.axis.parter.linear(tickdists
                                                                = [6 * 60,
                                                                   3 * 60]),
                              texter = time_of_day(base_time
                                                   = start_time_of_day_min),
                              painter = grid_painter),
        y = graph.axis.log(min = 1e-23,
                           max = 1e-17,
                           title = r"mass (kg)",
                           painter = grid_painter),
        key = graph.key.key(pos = "tr", hinside = 0,
                            symbolwidth = unit.v_cm))

    plot_data = [[] for s in aero_species]
    for [time, particle_list] in particle_history:
        particle = particle_by_id(particle_list, show_particles[i]["id"])
        if particle == None:
            continue
        for s in range(len(aero_species)):
            plot_data[s].append([time / 60,
                                 particle.mass(include = aero_species[s][1])])
    if max([len(d) for d in plot_data]) == 0:
        raise Exception("Particle ID not found: %d" % show_particles[i]["id"])

    for s in range(len(aero_species)):
        plot_data[s].sort()
        plot_data[s] = [[time, value] for [time, value] in plot_data[s]
                        if value > 0.0]
        if aero_species[s][0] == "":
            label = tex_species(aero_species[s][1][0])
        else:
            label = aero_species[s][0]
        if len(plot_data[s]) > 0:
            g.plot(graph.data.points(plot_data[s], x = 1, y = 2, title = label),
                   styles = [graph.style.line(lineattrs
                                              = [line_style_list[s],
                                                 style.linewidth.THick])])
        #label_plot_line(g, plot_data[s], 18 * 60.0,
        #            label, [0, 1], 1 * unit.v_mm)

    g.doaxes()
    g.dodata()
    boxed_text(g, 0.04, 0.9, show_particles[i]["box label"])
    g.writePDFfile("%s_%s.pdf" % (out_prefix, show_particles[i]["suffix"]))
