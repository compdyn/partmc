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
netcdf_pattern_wc = r"urban_plume_wc_state_0001_([0-9]{8})\.nc"

netcdf_dir_nc = "out"
netcdf_pattern_nc = r"urban_plume_nc_state_0001_([0-9]{8})\.nc"

max_val = 4.0

diameter_axis_min = 0.01
diameter_axis_max = 1.0
num_diameter_bins = 70

times_hour = {"g11": 1,
              "g12" : 5,
              "g21" : 7,
              "g22" : 24}

grid_v_space = 0.7
grid_h_space = 0.5
grid_graph_width = 6.3

show_particles = [
#    {"id": 106390, "suffix": "1", "label": "P1",
#     "label pos": [0.9, 0.4], "box label": "particle P1"},
#    {"id": 195377, "suffix": "2", "label": "P2",
#     "label pos": [0.1, 0.5], "box label": "particle P2"},
    {"id": 108139, "suffix": "1", "label": "P1",
     "label pos": [0.9, 0.4], "box label": "particle P1"},
    {"id": 192536, "suffix": "2", "label": "P2",
     "label pos": [0.1, 0.5], "box label": "particle P2"},
    {"id": 549273, "suffix": "3", "label": "P3",
     "label pos": [0.9, 0.7], "box label": "particle P3"},
    ]

def make_2x2_graph_grid(y_axis_label):
    c = canvas.canvas()
    g21 = c.insert(graph.graphxy(
        width = grid_graph_width,
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
                           title = r'dry diameter ($\mu$m)'),
        y = graph.axis.linear(min = 0,
                              max = 100,
                              title = y_axis_label,
                              texter = graph.axis.texter.decimal(suffix
                                                                 = r"\%"))))
    g11 = c.insert(graph.graphxy(
        width = grid_graph_width,
        ypos = g21.height + grid_v_space,
        x = graph.axis.linkedaxis(g21.axes["x"]),
        y = graph.axis.linear(min = 0,
                              max = 100,
                              title = y_axis_label,
                              texter = graph.axis.texter.decimal(suffix
                                                                 = r"\%"))))
    g22 = c.insert(graph.graphxy(
        width = grid_graph_width,
        xpos = g21.width + grid_h_space,
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
                           title = r'dry diameter ($\mu$m)'),
        y = graph.axis.linkedaxis(g21.axes["y"])))
    g12 = c.insert(graph.graphxy(
        width = grid_graph_width,
        xpos = g11.width + grid_h_space,
        ypos = g22.height + grid_v_space,
        x = graph.axis.linkedaxis(g22.axes["x"]),
        y = graph.axis.linkedaxis(g11.axes["y"])))
    return {"c": c,
            "g11" : g11,
            "g12" : g12,
            "g21" : g21,
            "g22" : g22}

def make_2x1_graph_grid(y_axis_label):
    c = canvas.canvas()
    
    g11 = c.insert(graph.graphxy(
        width = grid_graph_width,
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
                           title = r'dry diameter ($\mu$m)'),
        y = graph.axis.linear(min = 0,
                              max = 100,
                              title = y_axis_label,
                              texter = graph.axis.texter.decimal(suffix
                                                                 = r"\%"))))
    g21 = c.insert(graph.graphxy(
        width = grid_graph_width,
        xpos = g11.width + grid_h_space,
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
                           title = r'dry diameter ($\mu$m)'),
        y = graph.axis.linkedaxis(g11.axes["y"])))

    return {"c": c,
            "g11" : g11,
            "g21" : g21}

def write_time(g, env_state, extra_text = "", text_vpos = [0, 1],
               anchor_point_rel = [0, 1]):
    time_lst = time_of_day_string(env_state.start_time_of_day
                                  + env_state.elapsed_time)
    time_hour = int(env_state.elapsed_time / 3600.0)
    suffix = "s"
    if time_hour == 1:
        suffix = ""
    text = "%d hour%s (%s LST)%s" % (time_hour, suffix, time_lst, extra_text)
    (x_g, y_g) = g.vpos(text_vpos[0], text_vpos[1])
    boxed_text_g(g, text, x_g, y_g, anchor_point_rel = anchor_point_rel)

def write_time_outside(g, env_state, extra_text = ""):
    write_time(g, env_state, extra_text)

def write_time_inside(g, env_state, extra_text = ""):
    write_time(g, env_state, extra_text, anchor_point_rel = [0, 0])

def write_text_outside(g, text):
    (x_g, y_g) = g.vpos(0, 1)
    boxed_text_g(g, text, x_g, y_g, anchor_point_rel = [0, 1])

def write_text_inside(g, text):
    (x_g, y_g) = g.vpos(0, 1)
    boxed_text_g(g, text, x_g, y_g, anchor_point_rel = [0, 0])
