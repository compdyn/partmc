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

netcdf_dir_wc = "out"
netcdf_pattern_wc = r"urban_plume_state_0001_([0-9]{8})\.nc"

netcdf_dir_nc = "out"
netcdf_pattern_nc = r"urban_plume_state_0001_([0-9]{8})\.nc"

max_val = 4.0

diameter_axis_min = 2e-3
diameter_axis_max = 1e0

times_hour = {"g11": 1,
              "g12" : 6,
              "g21" : 9,
              "g22" : 24}

show_particles = [
    {"id": 1285, "suffix": "wet_diesel", "label": "wet diesel",
     "label pos": [0.4, 55], "box label": "wet diesel particle"},
    {"id": 2511, "suffix": "dry_diesel", "label": "dry diesel",
     "label pos": [0.006, 55], "box label": "dry diesel particle"},
    {"id": 5027, "suffix": "late_diesel", "label": "late diesel",
     "label pos": [0.4, 75], "box label": "late diesel particle"},
    ]

def make_4x4_graph_grid(y_axis_label):
    v_space = 0.5
    h_space = 0.5
    graph_width = 6.3
    
    graphs = {}
    c = canvas.canvas()
    
    g21 = c.insert(graph.graphxy(
        width = graph_width,
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
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
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
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

def corner_boxed_text(g, text):
    boxed_text(g, 0.04, 0.9, text)

def write_time(g, env_state):
    time_lst = time_of_day_string(env_state.start_time_of_day
                                  + env_state.elapsed_time)
    time_hour = int(env_state.elapsed_time / 3600.0)
    suffix = "s"
    if time_hour == 1:
        suffix = ""
    corner_boxed_text(g, "%d hour%s (%s LST)" % (time_hour, suffix, time_lst))
