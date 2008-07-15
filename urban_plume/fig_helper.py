#!/usr/bin/env python
# Copyright (C) 2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *

max_val = 4.0

times_hour = {"g11": 1,
              "g12" : 3,
              "g21" : 6,
              "g22" : 8}

show_particles = [[1049, [0.02, 75], "wet diesel"],
                  [1279, [0.02, 15], "wet gasoline"],
                  [2329, [0.7, 75], "dry diesel"],
                  [2315, [0.7, 15], "dry gasoline"],
                  ]
show_particles = []

def get_time_filename_list(dir, file_pattern):
    time_filename_list = []
    filenames = os.listdir(dir)
    if filenames == []:
        raise Exception("No files in %s match %s" % (dir, file_pattern))
    file_re = re.compile(file_pattern)
    for filename in filenames:
        match = file_re.search(filename)
        if match:
            output_key = match.group(1)
            netcdf_filename = os.path.join(dir, filename)
            ncf = NetCDFFile(netcdf_filename)
            env_state = env_state_t(ncf)
            time_filename_list.append([env_state.elapsed_time,
                                       netcdf_filename])
            ncf.close()
    time_filename_list.sort()
    return time_filename_list

def file_filename_at_time(time_filename_list, search_time):
    min_diff = abs(search_time - time_filename_list[0][0])
    min_filename = time_filename_list[0][1]
    for [time, filename] in time_filename_list[1:]:
        diff = abs(search_time - time)
        if diff < min_diff:
            min_diff = diff
            min_filename = filename
    return min_filename

def make_4x4_graph_grid(y_axis_label):
    v_space = 0.5
    h_space = 0.5
    graph_width = 6.3
    
    graphs = {}
    c = canvas.canvas()
    
    g21 = c.insert(graph.graphxy(
        width = graph_width,
        x = graph.axis.log(min = 2.e-3,
                           max = 1.e+0,
                           title = r'dry diameter ($\mu$m)'),
        y = graph.axis.linear(min = 0,
                              max = 100,
                              title = y_axis_label,
                              texter = graph.axis.texter.decimal(suffix
                                                                 = r"\%"))))
    g11 = c.insert(graph.graphxy(
        width = graph_width,
        ypos = g21.height + v_space,
        x = graph.axis.linkedaxis(g21.axes["x"]),
        y = graph.axis.linear(min = 0,
                              max = 100,
                              title = y_axis_label,
                              texter = graph.axis.texter.decimal(suffix
                                                                 = r"\%"))))
    g22 = c.insert(graph.graphxy(
        width = graph_width,
        xpos = g21.width + h_space,
        x = graph.axis.log(min = 2.e-3,
                           max = 1.e+0,
                           title = r'dry diameter ($\mu$m)'),
        y = graph.axis.linkedaxis(g21.axes["y"])))
    g12 = c.insert(graph.graphxy(
        width = graph_width,
        xpos = g11.width + h_space,
        ypos = g22.height + v_space,
        x = graph.axis.linkedaxis(g22.axes["x"]),
        y = graph.axis.linkedaxis(g11.axes["y"])))
    return {"c": c,
            "g11" : g11,
            "g12" : g12,
            "g21" : g21,
            "g22" : g22}

