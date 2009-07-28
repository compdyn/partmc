#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math, random
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append(".")
from fig_helper import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *

out_prefix = "figs_aging/aging_aero_2d_age"

ss_axis_min = 1e-2
ss_axis_max = 1e2
num_ss_bins = 40

bc_max_val = 4.0
ss_max_val = 100.0

const = load_constants("../src/constants.f90")

time_filename_list = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
for color in [True, False]:
    c = canvas.canvas()
    g1 = c.insert(graph.graphxy(
        width = grid_graph_width,
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
                           title = diameter_axis_label),
        y = graph.axis.log(min = ss_axis_min,
                              max = ss_axis_max,
                              title = r"critical supersaturation $S_{\rm c}$ (\%)")))
    g2 = c.insert(graph.graphxy(
            width = grid_graph_width,
            ypos = g1.height + grid_v_space,
            x = graph.axis.linkedaxis(g1.axes["x"]),
            y = graph.axis.linear(min = bc_axis_min,
                                  max = bc_axis_max,
                                  title = r"BC dry mass frac. $w_{{\rm BC},{\rm dry}}\ (\%)$",
                                  density = 1.2)))

    time = 24 * 3600.0
    filename = file_filename_at_time(time_filename_list, time)

    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()

    diameter = particles.dry_diameter() * 1e6
    critical_ss = (particles.kappa_rh(env_state, const) - 1.0) * 100.0
    comp_frac = particles.mass(include = ["BC"]) \
                / particles.mass(exclude = ["H2O"]) * 100
    age = (env_state.elapsed_time - particles.least_create_time) / 3600.0
    max_age = 24.0
    scaled_age = age / max_age

    point_size = 0.01

    ss_plot_data = zip(diameter, critical_ss, point_size * ones_like(diameter), scaled_age)
    ss_plot_data = [(x, y, r, clr) for (x, y, r, clr) in ss_plot_data
                    if x >= diameter_axis_min
                    and x <= diameter_axis_max
                    and y >= ss_axis_min
                    and y <= ss_axis_max]
    random.shuffle(ss_plot_data)

    bc_plot_data = zip(diameter, comp_frac, point_size * ones_like(diameter), scaled_age)
    bc_plot_data = [(x, y, r, clr) for (x, y, r, clr) in bc_plot_data
                    if x >= diameter_axis_min
                    and x <= diameter_axis_max
                    and y >= bc_axis_min
                    and y <= bc_axis_max]
    random.shuffle(bc_plot_data)

    if color:
        palette = rainbow_palette
    else:
        palette = nonlinear_gray_palette

    for g in [g1, g2]:
        g.dolayout()

        for axisname in ["x", "y"]:
            for t in g.axes[axisname].data.ticks:
                if t.ticklevel is not None:
                    g.stroke(g.axes[axisname].positioner.vgridpath(t.temp_v),
                             [style.linestyle.dotted])
                    
    for (plot_data, g) in [(ss_plot_data, g1),
                           (bc_plot_data, g2)]:
        for [x, y, r, clr] in plot_data:
            x_pt, y_pt = g.pos(x, y)
            g.draw(path.circle(x_pt, y_pt, r),
                   [deco.filled([palette.getcolor(clr)])])

    for g in [g1, g2]:
        write_time(g, env_state)

        g.dodata()
        g.doaxes()

        add_canvas_color_bar(
            g,
            min = 0.0,
            max = max_age,
            xpos = g.xpos + g.width + grid_h_space,
            ybottom = g.ypos,
            ytop = g.ypos + g.height,
            title = r"age (hours)",
            palette = palette,
            parter = graph.axis.parter.linear(tickdists = [6, 3]))


    if color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    c.writePDFfile(out_filename)
    if not color:
        print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
