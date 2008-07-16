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

aero_species = [
    {"species": ["NO3"], "plot": "g1", "label_time": 8, "label_pos": [1, 1]},
    {"species": ["NH4"], "plot": "g1", "label_time": 11, "label_pos": [1, 1]},
    {"species": ["OC"], "plot": "g1", "label_time": 10, "label_pos": [0, 0]},
    {"species": ["SO4"], "plot": "g2", "label_time": 9, "label_pos": [1, 1]},
    {"species": ["BC"], "plot": "g2", "label_time": 8, "label_pos": [0, 0]},
    {"species": ["ARO1", "ARO2", "ALK1", "OLE1"], "plot": "g2",
     "label": "SOA", "label_time": 8, "label_pos": [1, 1]},
    ]

out_filename = "figs/aero_time.pdf"

netcdf_dir = "out"
netcdf_pattern = r"urban_plume_state_0001_([0-9]{8})\.nc"

time_filename_list = get_time_filename_list(netcdf_dir, netcdf_pattern)
env_state = read_any(env_state_t, netcdf_dir, netcdf_pattern)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max([time for [time, filename] in time_filename_list]) / 60

c = canvas.canvas()

g2 = c.insert(graph.graphxy(
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
                          max = 10,
                          title = r"mass density ($\rm \mu g \, m^{-3}$)",
			  painter = grid_painter)))
g1 = c.insert(graph.graphxy(
    width = 6.7,
    ypos = g2.height + 0.5,
    x = graph.axis.linkedaxis(g2.axes["x"],
                              painter = graph.axis.painter.linked(gridattrs = [style.linestyle.dotted])),
    y = graph.axis.linear(min = 0.,
                          max = 30,
                          title = r"mass density ($\rm \mu g \, m^{-3}$)",
			  painter = grid_painter)))
plot_data = [[] for x in aero_species]
max_comp_vol = None
min_comp_vol = None
for [time, filename] in time_filename_list:
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    ncf.close()
    for i in range(len(aero_species)):
        masses = particles.mass(include = aero_species[i]["species"])
        mass_den = (masses / particles.comp_vol).sum()
        plot_data[i].append([time / 60.0, mass_den * 1e9])
    if max_comp_vol == None:
        max_comp_vol = particles.comp_vol.max()
    else:
        max_comp_vol = max(max_comp_vol, particles.comp_vol.max())
    if min_comp_vol == None:
        min_comp_vol = particles.comp_vol.min()
    else:
        min_comp_vol = min(min_comp_vol, particles.comp_vol.min())

print "max comp_vol = %g cm^3" % (max_comp_vol * 1e6)
print "min comp_vol = %g cm^3" % (min_comp_vol * 1e6)

graphs = {"g1": g1, "g2": g2}
line_counts = {"g1": 0, "g2": 0}
for i in range(len(aero_species)):
    graph_name = aero_species[i]["plot"]
    g = graphs[graph_name]
    g.plot(
        graph.data.points(plot_data[i], x = 1, y = 2),
        styles = [graph.style.line(lineattrs
                                   = [line_style_list[line_counts[graph_name]],
                                      style.linewidth.Thick])])
    line_counts[graph_name] += 1
    if "label" in aero_species[i].keys():
        label = aero_species[i]["label"]
    else:
        label = tex_species(aero_species[i]["species"][0])
    label_plot_line(g, plot_data[i], aero_species[i]["label_time"] * 60.0,
                    label, aero_species[i]["label_pos"], 1 * unit.v_mm)

c.writePDFfile(out_filename)
print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
