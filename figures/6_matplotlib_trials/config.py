#!/usr/bin/env python
# Copyright (C) 2008-2010 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math



netcdf_dir = os.path.join("../../scenarios/4_nucleate/out")
netcdf_indexed_patterns = [
    [1, r"^urban_plume_wc_0001_(00000001)\.nc$"],
    [2, r"^urban_plume_wc_0001_(00000007)\.nc$"],
    [3, r"^urban_plume_wc_0001_(00000013)\.nc$"],
    [4, r"^urban_plume_wc_0001_(00000019)\.nc$"],
    [5, r"^urban_plume_wc_0001_(00000025)\.nc$"],
    [6, r"^urban_plume_wc_0001_(00000031)\.nc$"],
    [7, r"^urban_plume_wc_0001_(00000037)\.nc$"],
    [8, r"^urban_plume_wc_0001_(00000043)\.nc$"],
    [9, r"^urban_plume_wc_0001_(00000049)\.nc$"],
    [10, r"^urban_plume_wc_0001_(00000055)\.nc$"],
    [11, r"^urban_plume_wc_0001_(00000061)\.nc$"],
    [12, r"^urban_plume_wc_0001_(00000067)\.nc$"],
    [13, r"^urban_plume_wc_0001_(00000073)\.nc$"],
    [14, r"^urban_plume_wc_0001_(00000079)\.nc$"],
    [15, r"^urban_plume_wc_0001_(00000085)\.nc$"],
    [16, r"^urban_plume_wc_0001_(00000091)\.nc$"],
    [17, r"^urban_plume_wc_0001_(00000097)\.nc$"],
    [18, r"^urban_plume_wc_0001_(00000103)\.nc$"],
    [19, r"^urban_plume_wc_0001_(00000109)\.nc$"],
    [20, r"^urban_plume_wc_0001_(00000115)\.nc$"],
    [21, r"^urban_plume_wc_0001_(00000121)\.nc$"],
    [22, r"^urban_plume_wc_0001_(00000127)\.nc$"],
    [23, r"^urban_plume_wc_0001_(00000133)\.nc$"],
    [24, r"^urban_plume_wc_0001_(00000139)\.nc$"],
    [25, r"^urban_plume_wc_0001_(00000145)\.nc$"],
#    [1, r"^cond_01_ref_0001_([0-9]{8})\.nc$"],
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
bc_axis_max = 100
num_bc_bins = 100

oc_axis_min = 0
oc_axis_max = 100
num_oc_bins = 100

soa_axis_min = 0
soa_axis_max = 30
num_soa_bins = 100

inorg_axis_min = 0
inorg_axis_max = 100
num_inorg_bins = 100
