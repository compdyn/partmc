#!/usr/bin/env python
# Copyright (C) 2008-2010 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math

netcdf_dir = os.path.join("../../scenarios/3_condense/out")
netcdf_indexed_patterns = [
    [1, r"^cond_01_ref_0001_([0-9]{8})\.nc$"],
#    [2, r"^cond_2_0001_([0-9]{8})\.nc$"],
#    [3, r"^cond_3_0001_([0-9]{8})\.nc$"],
#    [4, r"^cond_4_0001_([0-9]{8})\.nc$"],
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
