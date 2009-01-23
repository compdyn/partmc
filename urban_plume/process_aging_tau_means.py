#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, pyx
from numpy import *
sys.path.append("../tool")
from pmc_pyx import *
sys.path.append(".")
from process_aging_config import *

g = pyx.graph.graphxy(
    width = 20,
    x = pyx.graph.axis.linear(title = r"critical supersaturation",
                              painter = grid_painter),
    y = pyx.graph.axis.log(title = r"characteristic time (hours)",
                           painter = grid_painter),
    key = pyx.graph.key.key(pos = "tr"))

red = pyx.color.rgb(1, 0, 0)
magenta = pyx.color.rgb(1, 0, 1)
blue = pyx.color.rgb(0, 0, 1)
cyan = pyx.color.rgb(0, 1, 1)

solid = pyx.style.linestyle.solid
dashed = pyx.style.linestyle(pyx.style.linecap.round, pyx.style.dash([2, 2]))
dashdot = pyx.style.linestyle(pyx.style.linecap.round, pyx.style.dash([0, 2, 2, 2]))

data = [
    [tau_day_wc_mass, "tau day wc mass", red, dashdot],
    [tau_day_wc_cond_mass, "tau day wc cond mass", red, solid],
    [tau_day_nc_mass, "tau day nc mass", red, dashed],
    [tau_day_wc_num, "tau day wc num", magenta, dashdot],
    [tau_day_wc_cond_num, "tau day wc cond num", magenta, solid],
    [tau_day_nc_num, "tau day nc num", magenta, dashed],
    [tau_night_wc_mass, "tau night wc mass", blue, dashdot],
    [tau_night_wc_cond_mass, "tau night wc cond mass", blue, solid],
    [tau_night_nc_mass, "tau night nc mass", blue, dashed],
    [tau_night_wc_num, "tau night wc num", cyan, dashdot],
    [tau_night_wc_cond_num, "tau night wc cond num", cyan, solid],
    [tau_night_nc_num, "tau night nc num", cyan, dashed],
    ]

for [tau_data, label, color, pattern] in data:
    tau_data_hour = [t / 3600 for t in tau_data]
    g.plot(
        pyx.graph.data.points(zip(ss_active, tau_data_hour), x = 1, y = 2,
                              title = label),
        styles = [pyx.graph.style.line(lineattrs = [color, pattern])])

g.writePDFfile("%s/aging_tau.pdf" % data_prefix)

