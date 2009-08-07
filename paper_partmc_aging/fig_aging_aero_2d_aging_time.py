#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import random as py_random
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append(".")
from fig_helper import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *

out_prefix = "figs_aging/aging_aero_2d_aging_time"

const = load_constants("../src/constants.f90")

coag = True
bin = level_mid + 1

for plot_info in aging_time_infos:
    if coag:
        env_state = read_any(env_state_t, netcdf_dir_wc, netcdf_pattern_wc)
        coag_suffix = "wc"
    else:
        env_state = read_any(env_state_t, netcdf_dir_nc, netcdf_pattern_nc)
        coag_suffix = "nc"
    filename = os.path.join(aging_data_dir,
                            "particle_aging_%s_%s_plot_data_%08d.txt" % (plot_info["name"], coag_suffix, bin))
    value = loadtxt(filename)
    mask = where(value > 0.0, 1.0, 0.0)
    max_val = value.max()
    value_zero_to_max = where(value > 0.0, value, max_val)
    min_val = value_zero_to_max.min()
    log_value = where(value > 0.0, log(value), 0.0)
    value = log_value
    value = (value - log(min_val)) / (log(max_val) - log(min_val))
    value = value.clip(0.0, 1.0)

    plot_data = pmc_histogram_2d(value, diameter_axis, aging_time_axis, mask = mask)

    for color in [True, False]:
        g = graph.graphxy(
            width = grid_graph_width,
            x = graph.axis.log(min = diameter_axis_min,
                               max = diameter_axis_max,
                               title = r"dry diameter at emission $D\ (\rm\mu m)$"),
            y = graph.axis.linear(min = aging_time_axis.min,
                                  max = aging_time_axis.max,
                                  parter = graph.axis.parter.linear(tickdists = [6, 3]),
                                  title = r"aging time $t_{\rm aging}\ (\rm h)$"))

        if color:
            palette = rainbow_palette
        else:
            palette = gray_palette

        g.dolayout()

        for axisname in ["x", "y"]:
            for t in g.axes[axisname].data.ticks:
                if t.ticklevel is not None:
                    g.stroke(g.axes[axisname].positioner.vgridpath(t.temp_v),
                             [style.linestyle.dotted])

        g.plot(graph.data.points(plot_data,
                                 xmin = 1, xmax = 2, ymin = 3, ymax = 4,
                                 color = 5),
               styles = [hsb_rect(palette)])

        g.dodata()
        g.doaxes()

        write_text_outside(g, r"critical supersaturation $S_{\rm c} = %.1f\%%$" % (level_mid_value * 100))

        if plot_info["criterion"] == "none":
            plot_info_text = "all particles"
        elif plot_info["criterion"] == "emission_time":
            start_lst = time_of_day_string(env_state.start_time_of_day
                                           + plot_info["time_start"],
                                           separator = ":")
            end_lst = time_of_day_string(env_state.start_time_of_day
                                         + plot_info["time_end"],
                                         separator = ":")
            plot_info_text = "particles emitted between %s LST and %s LST" % (start_lst, end_lst)
        elif plot_info["criterion"] == "aging_time":
            start_lst = time_of_day_string(env_state.start_time_of_day
                                           + plot_info["time_start"],
                                           separator = ":")
            end_lst = time_of_day_string(env_state.start_time_of_day
                                         + plot_info["time_end"],
                                         separator = ":")
            plot_info_text = "particles that aged between %s LST and %s LST" % (start_lst, end_lst)
        (x_g, y_g) = g.vpos(0, 1.1)
        boxed_text_g(g, plot_info_text, x_g, y_g, anchor_point_rel = [0, 1])

        add_canvas_color_bar(
            g,
            log_scale = True,
            min = min_val,
            max = max_val,
            xpos = g.xpos + g.width + grid_h_space,
            ybottom = g.ypos,
            ytop = g.ypos + g.height,
            title = r"normalized number conc.",
            palette = palette)

        if color:
            out_filename = "%s_%s_%s_color.pdf" % (out_prefix, plot_info["name"], coag_suffix)
        else:
            out_filename = "%s_%s_%s_bw.pdf" % (out_prefix, plot_info["name"], coag_suffix)
        g.writePDFfile(out_filename)
        if not color:
            print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
            print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
