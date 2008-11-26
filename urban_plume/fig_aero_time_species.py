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

aero_species = [
    {"species": ["NO3"], "plot": "g1", "label_time": 8, "label_pos": [1, 1]},
    {"species": ["NH4"], "plot": "g1", "label_time": 8, "label_pos": [0, 0]},
    {"species": ["OC"], "plot": "g1", "label_time": 5, "label_pos": [1, 1]},
    {"species": ["SO4"], "plot": "g2", "label_time": 6, "label_pos": [1, 1]},
    {"species": ["BC"], "plot": "g2", "label_time": 10, "label_pos": [1, 1]},
    {"species": ["ARO1", "ARO2", "ALK1", "OLE1"], "plot": "g2",
     "label": "SOA", "label_time": 8, "label_pos": [0, 0]},
    ]

out_prefix = "figs/aero_time_species"

time_filename_list = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
start_time_of_day_min = env_state.start_time_of_day / 60
max_time_min = max([time for [time, filename, key] in time_filename_list]) / 60

for use_color in [True, False]:
    c = canvas.canvas()

    g2 = c.insert(graph.graphxy(
        width = 6.8,
        x = graph.axis.linear(min = 0.,
                              max = max_time_min,
                              parter = graph.axis.parter.linear(tickdists
                                                                = [6 * 60, 3 * 60]),
                              texter = time_of_day(base_time
                                                   = start_time_of_day_min),
                              title = "local standard time (LST) (hours:minutes)",
                              painter = grid_painter),
        y = graph.axis.linear(min = 0.0,
                              max = 8.0,
                              density = 1.5, # hack to increase number of ticks
                              title = r"mass conc. ($\rm \mu g \, m^{-3}$)",
                              painter = grid_painter)))
    g1 = c.insert(graph.graphxy(
        width = 6.7,
        ypos = g2.height + 0.5,
        x = graph.axis.linkedaxis(g2.axes["x"],
                                  painter = graph.axis.painter.linked(gridattrs = [style.linestyle.dotted])),
        y = graph.axis.linear(min = 0.0,
                              max = 30.0,
                              title = r"mass conc. ($\rm \mu g \, m^{-3}$)",
                              painter = grid_painter)))

    g1.doaxes()
    g2.doaxes()

    plot_data = [[] for x in aero_species]
    max_comp_vol = None
    min_comp_vol = None
    max_n_particles = None
    min_n_particles = None
    found_dry_diesel = False
    found_water_transition = False
    found_dry_diesel_with_nitrate = False
    for [time, filename, key] in time_filename_list:
        ncf = NetCDFFile(filename)
        particles = aero_particle_array_t(ncf)
        env_state = env_state_t(ncf)
        ncf.close()

        bc_oc_frac = particles.mass(include = ["BC"]) \
                     / particles.mass(include = ["BC", "OC"]) * 100
        water_frac = particles.mass(include = ["H2O"]) \
                     / particles.mass() * 100
        nitrate_frac = particles.mass(include = ["NO3"]) \
                       / particles.mass(exclude = ["H2O"]) * 100
        soot_water_frac = water_frac[bc_oc_frac > 2]
        wet_soot_water_frac = soot_water_frac[soot_water_frac > 0.0]
        n_orig_part = array([int(n) for n in particles.n_orig_part])
        if len(soot_water_frac) > 0:
            fraction_wet_soot = len(wet_soot_water_frac) \
                                / float(len(soot_water_frac)) * 100
            if not found_water_transition:
                if (time > 6 * 3600) and (fraction_wet_soot > 50):
                    found_water_transition = True
                    if not use_color:
                        print ("Water transition after %g seconds (at %s LST)"
                               " with RH = %g%%") \
                               % (env_state.elapsed_time,
                                  time_of_day_string(env_state.elapsed_time
                                                     + env_state.start_time_of_day),
                                  env_state.relative_humidity * 100)
        if not found_dry_diesel_with_nitrate:
            if ((bc_oc_frac > 60) & (water_frac == 0.0) \
                & (nitrate_frac > 0.0) & (n_orig_part == 1)).any():
                found_dry_diesel_with_nitrate = True
                if not use_color:
                    print ("First dry diesel particle with nitrate after %g seconds"
                           " (at %s LST) with RH = %g%%") \
                           % (env_state.elapsed_time,
                              time_of_day_string(env_state.elapsed_time
                                                 + env_state.start_time_of_day),
                              env_state.relative_humidity * 100)
        if not found_dry_diesel:
            if ((bc_oc_frac > 60) & (water_frac == 0.0)).any():
                found_dry_diesel = True
                if not use_color:
                    print ("First dry diesel particle after %g seconds (at %s LST)"
                           " with RH = %g%%") \
                           % (env_state.elapsed_time,
                              time_of_day_string(env_state.elapsed_time
                                                 + env_state.start_time_of_day),
                              env_state.relative_humidity * 100)

        for i in range(len(aero_species)):
            masses = particles.mass(include = aero_species[i]["species"])
            mass_den = (masses / particles.comp_vol).sum()
            plot_data[i].append([time / 60.0, mass_den * 1e9])
        if max_comp_vol == None:
            max_comp_vol = array(particles.comp_vol).max()
        else:
            max_comp_vol = max(max_comp_vol, array(particles.comp_vol).max())
        if min_comp_vol == None:
            min_comp_vol = array(particles.comp_vol).min()
        else:
            min_comp_vol = min(min_comp_vol, array(particles.comp_vol).min())
        if max_n_particles == None:
            max_n_particles = particles.n_particles
        else:
            max_n_particles = max(max_n_particles, particles.n_particles)
        if min_n_particles == None:
            min_n_particles = particles.n_particles
        else:
            min_n_particles = min(min_n_particles, particles.n_particles)

    if not use_color:
        if not found_water_transition:
            print "ERROR: did not find water transition"
        if not found_dry_diesel_with_nitrate:
            print "ERROR: did not find dry diesel with nitrate"
        if not found_dry_diesel:
            print "ERROR: did not find dry diesel"
        print "max comp_vol = %g cm^3" % (max_comp_vol * 1e6)
        print "min comp_vol = %g cm^3" % (min_comp_vol * 1e6)
        print "max n_particles = %d" % max_n_particles
        print "min n_particles = %d" % min_n_particles

    graphs = {"g1": g1, "g2": g2}
    line_counts = {"g1": 0, "g2": 0}
    for i in range(len(aero_species)):
        graph_name = aero_species[i]["plot"]
        g = graphs[graph_name]
        if use_color:
            style_attr = color_list[line_counts[graph_name]]
        else:
            style_attr = line_style_list[line_counts[graph_name]]
        g.plot(
            graph.data.points(plot_data[i], x = 1, y = 2),
            styles = [graph.style.line(lineattrs
                                       = [style_attr,
                                          style.linewidth.Thick])])
        line_counts[graph_name] += 1
        if "label" in aero_species[i].keys():
            label = aero_species[i]["label"]
            print_label = aero_species[i]["label"]
        else:
            label = tex_species(aero_species[i]["species"][0])
            print_label = aero_species[i]["species"][0]
        label_plot_line_boxed(g, plot_data[i], aero_species[i]["label_time"] * 60.0,
                        label, aero_species[i]["label_pos"])

        a_data = array([v for [t,v] in plot_data[i]])
        j = a_data.argmax()
        if not use_color:
            print "%s max = %g ug/m^3 at %s LST, initial = %g ug/m^3" \
                  % (print_label, plot_data[i][j][1],
                     time_of_day_string(plot_data[i][j][0] * 60
                                        + env_state.start_time_of_day),
                     plot_data[i][0][1])

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    c.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
