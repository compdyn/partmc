#!/usr/bin/env python
# Copyright (C) 2007-2010 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import Scientific.IO.NetCDF
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rc('text', usetex = True)
matplotlib.rc('xtick.major', pad = 8)
matplotlib.rc('ytick.major', pad = 8)
matplotlib.rc('xtick', labelsize = 10)
matplotlib.rc('legend', fontsize = 10, borderpad = 0.7, borderaxespad = 1)
matplotlib.rc('font', size = 10, family = "serif", serif = ["Computer Modern Roman"])
matplotlib.rc('lines', linewidth = 0.5)
matplotlib.rc('patch', linewidth = 0.5)
matplotlib.rc('axes', linewidth = 0.5)

out_filename_max_ss = "figs/mp_time_max_ss.pdf"
out_filename_active = "figs/mp_time_active.pdf"

particle_data = [
    [1,  0.1741616308,  41.6524],
    [7,  0.09059989756, 39.5619],
    [15, 0.07041647553, 43.5748],
    [24, 0.09695690379, 56.1   ],
    ]

size_avg_data = [
    [1,  0.1714984925,  40.7481],
    [7,  0.08763639475, 36.49  ],
    [15, 0.06879674967, 37.5584],
    [24, 0.09493217247, 57.6948],
    ]

comp_avg_data = [
    [1,  0.1422643991,  50.6337],
    [7,  0.07910420296, 45.3604],
    [15, 0.06421621256, 51.6355],
    [24, 0.08472553843, 65.4554],
    ]

golden_ratio = (1 + math.sqrt(5)) / 2
ax_width = 8 / 2.54 # inches
ax_height = ax_width / golden_ratio
ax_left_margin   = 0.7 # inches
ax_right_margin  = 0.2 # inches
ax_bottom_margin = 0.5 # inches
ax_top_margin    = 0.2 # inches
fig_width = ax_left_margin + ax_width + ax_right_margin
fig_height = ax_bottom_margin + ax_height + ax_top_margin
ax_left_margin_fraction = ax_left_margin / fig_width
ax_width_fraction = ax_width / fig_width
ax_bottom_margin_fraction = ax_bottom_margin / fig_height
ax_height_fraction = ax_height / fig_height

fig_max_ss = plt.figure()
#ax_max_ss = fig_max_ss.add_subplot(111)
ax_max_ss = fig_max_ss.add_axes([ax_left_margin_fraction,
                                 ax_bottom_margin_fraction,
                                 ax_width_fraction,
                                 ax_height_fraction])

fig_active = plt.figure()
ax_active = fig_active.add_subplot(111)

max_ss_lines = []
active_lines = []
labels = []

for (i_data, (data, name)) in enumerate(zip([particle_data, size_avg_data, comp_avg_data],
                                            [r"particle", r"size-avg", r"comp-avg"])):
    x = [d[0] for d in data]
    y = [d[1] for d in data]
    max_ss_lines.append(ax_max_ss.plot(x, y))
    y = [d[2] for d in data]
    active_lines.append(ax_active.plot(x, y))
    labels.append(name)

ax_max_ss.legend(max_ss_lines, labels, 'upper right')
ax_active.legend(active_lines, labels, 'upper left')

fig_max_ss.set_figwidth(fig_width)
fig_max_ss.set_figheight(fig_height)
ax_max_ss.grid(True)
ax_max_ss.grid(True, which = 'minor')
ax_max_ss.minorticks_on()

ax_max_ss.set_xticks([0, 6, 12, 18, 24])
ax_max_ss.set_xticks([3, 9, 15, 21], minor = True)
ax_max_ss.set_xbound(0, 24)
ax_max_ss.set_xlabel(r"time (h)")
ax_max_ss_x = ax_max_ss.get_xaxis()
ax_max_ss_y = ax_max_ss.get_yaxis()
#ax_max_ss_x.set_label_coords(0.5, -0.1)
ax_max_ss_x.labelpad = 8
ax_max_ss_y.labelpad = 8

ax_max_ss.set_ybound(0, 0.2)
ax_max_ss.set_yticks([0.025, 0.075, 0.125, 0.175], minor = True)

ax_max_ss.set_ylabel(r"maximum supersaturation (\%)")

#fig_max_ss.subplots_adjust(left = 0.2, bottom = 0.3)

fig_active.set_figwidth(4)

fig_max_ss.savefig(out_filename_max_ss)
fig_active.savefig(out_filename_active)
