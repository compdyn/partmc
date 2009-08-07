#!/usr/bin/env python
# Copyright (C) 2008, 2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *

text.set(mode="latex")
#text.set(fontmaps="download35.map") # embed times in the pdfs
#text.set(mode="latex",usefiles=["spam.aux"],texdebug="spam.debug")
from pmc_pyx import *

text.preamble(r"""\usepackage{times}
%\usepackage{mathptmx} % symbol math
\renewcommand{\normalsize}{\fontsize{9}{11}\selectfont}""")

#text.preamble(r"""\renewcommand{\sfdefault}{phv}
#\renewcommand{\familydefault}{\sfdefault}
#\renewcommand{\normalsize}{\fontsize{9}{11}\selectfont}""")

netcdf_dir_wc = "../urban_plume/out"
netcdf_pattern_wc = r"^urban_plume_wc_state_0001_([0-9]{8})\.nc$"

netcdf_dir_nc = "../urban_plume/out"
netcdf_pattern_nc = r"^urban_plume_nc_state_0001_([0-9]{8})\.nc$"

aging_data_dir = "out"
n_level_bin = 100
ss_active_axis = pmc_log_axis(0.0001, 0.1, n_level_bin)
level_low = int(ss_active_axis.closest_edge(array(0.001)))
level_mid = int(ss_active_axis.closest_edge(array(0.003)))
level_high = int(ss_active_axis.closest_edge(array(0.006)))
level_low_value = ss_active_axis.edge(level_low)
level_mid_value = ss_active_axis.edge(level_mid)
level_high_value = ss_active_axis.edge(level_high)

movie_stroke_width = 0.01

smooth_window_len = 60
grey_level = 0.2
grey_level_dark = 0.2
grey_level_light = 0.1

max_val = 4.0

diameter_axis_min = 0.01
diameter_axis_max = 1.0
num_diameter_bins = 70
diameter_axis_label = r'dry diameter $D\ (\rm\mu m)$'
diameter_axis = pmc_log_axis(min = diameter_axis_min, max = diameter_axis_max,
                             n_bin = num_diameter_bins)
diameter_axis_hi = pmc_log_axis(min = diameter_axis.min,
                                max = diameter_axis.max,
                                n_bin = diameter_axis.n_bin * 2)

bc_axis_min = 0
bc_axis_max = 80
num_bc_bins = 40

aging_time_axis = pmc_linear_axis(min = 0,
                                  max = 18,
                                  n_bin = 36)
aging_time_axis_hi = pmc_linear_axis(min = aging_time_axis.min,
                                     max = aging_time_axis.max,
                                     n_bin = aging_time_axis.n_bin * 2)
aging_time_axis_label = r'aging time $t_{\rm age}\rm\ (h)$'

aging_time_infos = [
    {"name": "all",
     "criterion": "none", "time_start" = -1, "time_end" = -1,
     "title": "all particles",
     },
    {"name": "emit_morning",
     "criterion": "emission_time", "time_start" = 5 * 3600.0, "time_end" = 6 * 3600.0,
     "title": None,
     },
    {"name": "emit_evening",
     "criterion": "emission_time", "time_start" = 11 * 3600.0, "time_end" = 12 * 3600.0,
     "title": None,
     },
    {"name": "age_day",
     "criterion": "aging_time", "time_start" = 6 * 3600.0, "time_end" = 9 * 3600.0,
     "title": None,
     },
    {"name": "age_night",
     "criterion": "aging_time", "time_start" = 12 * 3600.0, "time_end" = 22 * 3600.0,
     "title": None,
     },
    ]

times_hour = {"g11": 1,
              "g12" : 5,
              "g21" : 7,
              "g22" : 24}

grid_v_space = 0.7
grid_h_space = 0.5
grid_graph_width = 6.45

show_particles = [
    {"id": 467312, "suffix": "1", "label": "P1",
     "label pos": [0.9, 0.4], "box label": "particle P1"},
    {"id": 467197, "suffix": "2", "label": "P2",
     "label pos": [0.1, 0.5], "box label": "particle P2"},
    ]

def make_2x2_graph_grid(y_axis_label, y_min = bc_axis_min, y_max = bc_axis_max,
                        with_y_percent = False, with_x_percent = False,
                        y_log = False, x_log = True, with_key = False,
                        x_axis_label = diameter_axis_label,
                        x_min = diameter_axis_min, x_max = diameter_axis_max,
                        y_density = 1.2):
    c = canvas.canvas()
    if with_y_percent:
        y_texter = graph.axis.texter.decimal(suffix = r"\%")
    else:
        y_texter = graph.axis.texter.mixed()
    if with_x_percent:
        x_texter = graph.axis.texter.decimal(suffix = r"\%")
    else:
        x_texter = graph.axis.texter.mixed()
    if y_log:
        y = graph.axis.log(min = y_min,
                           max = y_max,
                           title = y_axis_label,
                           texter = y_texter,
                           density = y_density)
    else:
        y = graph.axis.linear(min = y_min,
                              max = y_max,
                              title = y_axis_label,
                              texter = y_texter,
                              density = y_density)
    if x_log:
        x = graph.axis.log(min = x_min,
                           max = x_max,
                           title = x_axis_label,
                           texter = x_texter)
    else:
        x = graph.axis.linear(min = x_min,
                              max = x_max,
                              title = x_axis_label,
                              texter = x_texter)
    g21 = c.insert(graph.graphxy(
        width = grid_graph_width,
        x = x,
        y = y))
    if with_key:
        key = graph.key.key(pos = "tr", vinside = 0, columns = 1)
    else:
        key = None
    g11 = c.insert(graph.graphxy(
            width = grid_graph_width,
            ypos = g21.height + grid_v_space,
            x = graph.axis.linkedaxis(g21.axes["x"]),
            y = y,
            key = key))
    g22 = c.insert(graph.graphxy(
        width = grid_graph_width,
        xpos = g21.width + grid_h_space,
        x = x,
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

def make_2x1_graph_grid(y_axis_label, y_density = 1.2):
    c = canvas.canvas()
    
    g11 = c.insert(graph.graphxy(
        width = grid_graph_width,
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
                           title = diameter_axis_label),
        y = graph.axis.linear(min = bc_axis_min,
                              max = bc_axis_max,
                              title = y_axis_label,
                              density = y_density)))
    g21 = c.insert(graph.graphxy(
        width = grid_graph_width,
        xpos = g11.width + grid_h_space,
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
                           title = diameter_axis_label),
        y = graph.axis.linkedaxis(g11.axes["y"])))

    return {"c": c,
            "g11" : g11,
            "g21" : g21}

def make_1x1_graph_grid(y_axis_label):
    g = graph.graphxy(
        width = grid_graph_width,
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
                           title = diameter_axis_label),
        y = graph.axis.linear(min = bc_axis_min,
                              max = bc_axis_max,
                              density = 1.2,
                              title = y_axis_label))
    return g

def write_time(g, env_state, extra_text = "", text_vpos = [0, 1],
               anchor_point_rel = [0, 1],
               with_hours = True):
    time_lst = time_of_day_string(env_state.start_time_of_day
                                  + env_state.elapsed_time,
                                  separator = ":")
    time_hour = int(env_state.elapsed_time / 3600.0)
    suffix = "s"
    if time_hour == 1:
        suffix = ""
    if with_hours:
        text = "%d hour%s (%s LST)%s" \
               % (time_hour, suffix, time_lst, extra_text)
    else:
        text = "%s LST%s" % (time_lst, extra_text)
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
