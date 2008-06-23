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

particles = [[1029, "wet_diesel"],
             [1892, "wet_gasoline"],
             [4116, "dry_diesel"],
             [4186, "dry_gasoline"],
             ]

netcdf_dir = "out"
netcdf_re = re.compile(r"urban_plume_0.5_3am_state_0001_[0-9]{8}\.nc")

aero_species_1 = [["", ["NO3"]],
                  ["", ["NH4"]],
                  ["", ["OC"]],
                  ["", ["H2O"]],
                  ["", ["SO4"]],
                  ["", ["BC"]],
                  ["SOA", ["ARO1", "ARO2", "ALK1", "OLE1"]],
                  ]

def interp(data, x):
    i = 0
    while data[i][0] < x:
        i += 1
    i -= 1
    if i < 0:
        raise Exception("data not in range: %s" % x)
    y = (x - data[i][0]) / (data[i+1][0] - data[i][0]) \
        * (data[i+1][1] - data[i][1]) + data[i][1]
    return y

def make_fig(particle_id, description):
    data_1 = [[] for s in aero_species_1]
    data_dry_diameter = []
    data_wet_diameter = []
    filenames = os.listdir(netcdf_dir)
    for filename in filenames:
        if netcdf_re.search(filename):
            netcdf_filename = os.path.join(netcdf_dir, filename)
            print netcdf_filename
            ncf = NetCDFFile(netcdf_filename)
            particles = read_particles(ncf, ids = [particle_id])
            if len(particles) > 0:
                particle = particles[0]
                time = float(ncf.variables["time"].getValue()) / 60.0
                for i in range(len(aero_species_1)):
                    data_1[i].append([time,
                                      particle.mass(include
                                                    = aero_species_1[i][1])])
                data_dry_diameter.append([time, particle.dry_diameter() * 1e6])
                data_wet_diameter.append([time, particle.diameter() * 1e6])
            ncf.close()
    if len(data_1[0]) == 0:
        raise Exception("particle ID not found: %d" % particle_id)

    first_time = min([x[0] for x in data_1[0]]) - 1
    first_time_lst = first_time + 3 * 60
    first_time_string = "%02d:%02d" % (int((first_time_lst / 60) % 24),
                                       int(first_time_lst % 60))
    data_1_shifted = [[] for x in data_1]
    for i in range(len(aero_species_1)):
        data_1_shifted[i] = [[x[0] - first_time, x[1]] for x in data_1[i]]

    c = canvas.canvas()

    g2 = c.insert(graph.graphxy(
        width = 8,
        x = graph.axis.linear(min = 0,
                              max = 24 * 60,
                              title = r'time (LST)',
                              parter = graph.axis.parter.linear(tickdists
                                                                = [6 * 60,
                                                                   3 * 60]),
                              texter = time_of_day(base_time = 3 * 60),
                              painter = grid_painter),
        y = graph.axis.log(min = 1e-28,
                           max = 1e-16,
                           title = r"mass (kg)",
                           painter = grid_painter),
        y2 = graph.axis.log(title = r"diameter ($\rm \mu m$)"),
        y3 = graph.axis.linear(title = "relative humidity (1)",
                               texter = graph.axis.texter.decimal(suffix
                                                                  = r"\%")),
        y4 = graph.axis.linear(title = "temperature (K)"),
        key = graph.key.key(pos = "br", hinside = 0, hdist = 5 * unit.v_cm)))

    g1 = c.insert(graph.graphxy(
        width = 8,
        ypos = g2.height + 2,
        x = graph.axis.log(min = 1,
                           max = 2e3,
                           title = r'time since %s LST (minutes)'
                           % first_time_string,
                           painter = grid_painter),
        y = graph.axis.log(min = 1e-28,
                           max = 1e-16,
                           title = r"mass (kg)",
                           painter = grid_painter),
        y2 = graph.axis.log(title = r"diameter ($\rm \mu m$)"),
        y3 = graph.axis.linear(title = "relative humidity (1)",
                               texter = graph.axis.texter.decimal(suffix
                                                                  = r"\%")),
        y4 = graph.axis.linear(title = "temperature (K)"),
        key = graph.key.key(pos = "br", hinside = 0, hdist = 5 * unit.v_cm)))

    for i in range(len(aero_species_1)):
        data_1[i].sort()
        data_1[i] = [x for x in data_1[i] if x[1] > 0.0]
        data_1_shifted[i].sort()
        data_1_shifted[i] = [x for x in data_1_shifted[i] if x[1] > 0.0]
        if aero_species_1[i][0] == "":
            title = tex_species(aero_species_1[i][1][0])
        else:
            title = aero_species_1[i][0]
        g1.plot(graph.data.points(data_1_shifted[i], x = 1, y = 2,
                                  title = title),
               styles = [graph.style.line(lineattrs = [color_list[i]])])
        g2.plot(graph.data.points(data_1[i], x = 1, y = 2, title = title),
               styles = [graph.style.line(lineattrs = [color_list[i]])])
    data_dry_diameter.sort()
    g2.plot(graph.data.points(data_dry_diameter, x = 1, y2 = 2,
                              title = "dry diameter"),
            styles = [graph.style.line(lineattrs = [color.rgb.black,
                                                    style.linestyle.solid,
                                                    style.linewidth.THick])])
    data_wet_diameter.sort()
    g2.plot(graph.data.points(data_wet_diameter, x = 1, y2 = 2,
                              title = "wet diameter"),
            styles = [graph.style.line(lineattrs = [color.rgb.black,
                                                    style.linestyle.solid])])

    data_dry_diameter_shifted = [[x[0] - first_time, x[1]]
                                 for x in data_dry_diameter
                                 if x[0] > first_time]
    data_wet_diameter_shifted = [[x[0] - first_time, x[1]]
                                 for x in data_wet_diameter
                                 if x[0] > first_time]
    g1.plot(graph.data.points(data_dry_diameter_shifted, x = 1, y2 = 2,
                              title = "dry diameter"),
    styles = [graph.style.line(lineattrs = [color.rgb.black,
                                            style.linestyle.solid,
                                            style.linewidth.THick])])
    g1.plot(graph.data.points(data_wet_diameter_shifted, x = 1, y2 = 2,
                              title = "wet diameter"),
    styles = [graph.style.line(lineattrs = [color.rgb.black,
                                            style.linestyle.solid])])

    ncf = NetCDFFile("out/urban_plume_0.5_3am_0001.nc")
    data = pmc_var(ncf, "env_state", [])
    data.scale_dim("time", 1.0/60)
    
    temp_data = module_copy.deepcopy(data)
    temp_data.reduce([select("env", "temp")])
    rh_data = module_copy.deepcopy(data)
    rh_data.reduce([select("env", "rel_humid")])
    rh_data.scale(100.0)
    height_data = module_copy.deepcopy(data)
    height_data.reduce([select("env", "height")])

    g2.plot(graph.data.points(rh_data.data_center_list(),
                              x = 1, y3 = 2,
                              title = "relative humidity"),
            styles = [graph.style.line(lineattrs = [color.rgb.black,
                                                    style.linestyle.dashed])])
    g2.plot(graph.data.points(temp_data.data_center_list(),
                              x = 1, y4 = 2,
                              title = "temperature"),
            styles = [graph.style.line(lineattrs = [color.rgb.black,
                                                    style.linestyle.dashdotted])])
    rh_shifted = [[x[0] - first_time, x[1]] \
                 for x in rh_data.data_center_list() \
                 if x[0] > first_time]
    rh_shifted[0:0] = [[1, interp(rh_data.data_center_list(),
                                  first_time + 1)]]
    temp_shifted = [[x[0] - first_time, x[1]] \
                   for x in temp_data.data_center_list() \
                   if x[0] > first_time]
    temp_shifted[0:0] = [[1, interp(temp_data.data_center_list(),
                                    first_time + 1)]]
    g1.plot(graph.data.points(rh_shifted,
                              x = 1, y3 = 2,
                              title = "relative humidity"),
            styles = [graph.style.line(lineattrs = [color.rgb.black,
                                                    style.linestyle.dashed])])
    g1.plot(graph.data.points(temp_shifted,
                              x = 1, y4 = 2,
                              title = "temperature"),
            styles = [graph.style.line(lineattrs = [color.rgb.black,
                                                    style.linestyle.dashdotted])])

    ncf.close()

    print "RH at 06:36 =", interp(rh_data.data_center_list(), \
          216)
    
    c.writePDFfile("out/particle_%s.pdf" % description)

for p in particles:
    make_fig(p[0], p[1])
