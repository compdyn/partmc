#!/usr/bin/env python
# Copyright (C) 2008, 2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../../tool")
from pmc_data_nc import *

text.set(mode="latex")
#text.set(mode="latex",usefiles=["spam.aux"],texdebug="spam.debug")
from pmc_pyx import *

#text.preamble(r"""\usepackage{times}
#\renewcommand{\normalsize}{\fontsize{9}{11}\selectfont}""")

#text.preamble(r"""\renewcommand{\sfdefault}{phv}
#\renewcommand{\familydefault}{\sfdefault}
#\renewcommand{\normalsize}{\fontsize{9}{11}\selectfont}""")

netcdf_dir = os.path.join("../../new_cond/out")
netcdf_indexed_patterns = [
    [1, r"^cond_1_0001_([0-9]{8})\.nc$"],
    [2, r"^cond_2_0001_([0-9]{8})\.nc$"],
    [3, r"^cond_3_0001_([0-9]{8})\.nc$"],
    [4, r"^cond_4_0001_([0-9]{8})\.nc$"],
    ]

graph_width = 8
color_bar_offset = 0.5

diameter_axis_min = 0.01
diameter_axis_max = 1.0
num_diameter_bins = 70
diameter_axis_label = r'dry diameter $D\ (\rm\mu m)$'

bc_axis_min = 0
bc_axis_max = 80
num_bc_bins = 40

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
