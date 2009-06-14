#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
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
sys.path.append(".")
from fig_helper import *

out_prefix = "figs/aero_particles"

aero_species = [
    {"label": "", "species": ["BC"],
     "style": style.linestyle.dashed, "thickness": style.linewidth.thick,
     "color": color_list[0]},
    {"label": "", "species": ["OC"],
     "style": style.linestyle.dashed, "thickness": style.linewidth.THick,
     "color": color_list[1]},
    {"label": "", "species": ["NO3"],
     "style": style.linestyle.solid, "thickness": style.linewidth.THick,
     "color": color_list[2]},
    {"label": "", "species": ["NH4"],
     "style": style.linestyle.dotted, "thickness": style.linewidth.THick,
     "color": color_list[3]},
    {"label": "", "species": ["SO4"],
     "style": style.linestyle.dotted, "thickness": style.linewidth.thick,
     "color": color_list[4]},
    {"label": "SOA", "species": ["ARO1", "ARO2", "ALK1", "OLE1"],
     "style": style.linestyle.dashdotted, "thickness": style.linewidth.thick,
     "color": color_list[5]},
    {"label": "", "species": ["H2O"],
     "style": style.linestyle.solid, "thickness": style.linewidth.thick,
     "color": color.gray.black},
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

for use_color in [True, False]:
    c = canvas.canvas()

    graphs = {}

    graphs[0] = c.insert(graph.graphxy(
        width = 6.4,
        x = graph.axis.linear(min = 0,
                              max = max_time_min,
                              title = r'local standard time (LST) (hours and minutes)',
                              parter = graph.axis.parter.linear(tickdists
                                                                = [6 * 60,
                                                                   3 * 60]),
                              texter = time_of_day(base_time
                                                   = start_time_of_day_min,
                                                   separator = ""),
                              painter = grid_painter),
        y = graph.axis.log(min = 1e-14,
                           max = 1e-8,
                           title = r"mass ($\rm\mu g$)",
                           painter = grid_painter)))

    for i in range(1, len(show_particles)):
        if i == len(show_particles) - 1:
            key = graph.key.key(pos = "tc", vinside = 0, columns = 4)
            #symbolwidth = unit.v_cm)
        else:
            key = None
        graphs[i] = c.insert(graph.graphxy(
            width = 6.4,
            ypos = graphs[i-1].ypos + graphs[i-1].height + 0.5,
            x = graph.axis.linkedaxis(graphs[i-1].axes["x"],
                                      painter = linked_grid_painter),
            y = graph.axis.log(min = 1e-14,
                               max = 1e-8,
                               title = r"mass ($\rm\mu g$)",
                               painter = grid_painter),
            key = key))

    for i in range(len(show_particles)):
        g = graphs[len(show_particles) - i - 1]

        plot_data = [[] for s in aero_species]
        for [time, particle_list] in particle_history:
            particle = particle_by_id(particle_list, show_particles[i]["id"])
            if particle == None:
                continue
            for s in range(len(aero_species)):
                spec_mass = particle.mass(include = aero_species[s]["species"]) * 1e9 # kg to ug
                plot_data[s].append([time / 60, spec_mass])
        if max([len(d) for d in plot_data]) == 0:
            raise Exception("Particle ID not found: %d" % show_particles[i]["id"])

        for s in range(len(aero_species)):
            plot_data[s].sort()
            plot_data[s] = [[time, value] for [time, value] in plot_data[s]
                            if value > 0.0]
            if aero_species[s]["label"] == "":
                label = tex_species(aero_species[s]["species"][0])
            else:
                label = aero_species[s]["label"]
            if len(plot_data[s]) > 0:
                if use_color:
                    attrs = [aero_species[s]["color"],
                             style.linewidth.thick]
                else:
                    attrs = [aero_species[s]["style"],
                             aero_species[s]["thickness"]]
                g.plot(graph.data.points(plot_data[s], x = 1, y = 2, title = label),
                       styles = [graph.style.line(lineattrs = attrs)])

        min_time_min = min([plot_data[s][0][0] for s in range(len(aero_species))])
        if not use_color:
            print "%s emitted at %s LST" \
                  % (show_particles[i]["label"],
                     time_of_day_string(min_time_min * 60
                                        + env_state.start_time_of_day))

        g.doaxes()
        g.dodata()

        write_text_inside(g, show_particles[i]["box label"])

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    c.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
