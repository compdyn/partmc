#!/usr/bin/env python
# Copyright (C) 2007, 2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
sys.path.append(".")
from fig_helper import *

time_hour = 24

diam_num_min_max = [1e1, 1e5]
diam_mass_min_max = [1e-4, 1e4]
bc_num_min_max = [1e3, 1e7]
bc_mass_min_max = [1e0, 1e4]

out_prefix = "figs/aero_2d_proj"

diam_axis_label = r'dry diameter $D$ ($\rm\mu m$)'
bc_axis_label = r"BC dry mass frac. $w_{{\rm BC},{\rm dry}}$ ($\%$)"
diam_num_axis_label = r"number conc. $n(D)$ ($\rm cm^{-3}$)"
diam_mass_axis_label = r"mass conc. $m(D)$ ($\rm \mu g\,m^{-3}$)"
bc_num_axis_label = r"number conc. $n_{\rm BC,dry}(w)$ ($\rm cm^{-3}$)"
bc_mass_axis_label = r"mass conc. $m_{\rm BC,dry}(w)$ ($\rm \mu g\,m^{-3}$)"

def get_plot_data(filename):
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()

    diameter = particles.dry_diameter() * 1e6
    comp_frac = particles.mass(include = ["BC"]) \
                / particles.mass(exclude = ["H2O"]) * 100
    mass = particles.mass() * 1e9

    diam_axis = pmc_log_axis(min = diameter_axis_min, max = diameter_axis_max,
                             n_bin = num_diameter_bins)
    bc_axis = pmc_linear_axis(min = bc_axis_min, max = bc_axis_max,
                              n_bin = num_bc_bins)
    diam_bin = diam_axis.find(diameter)
    # hack to avoid landing just around the integer boundaries
    comp_frac *= (1.0 + 1e-12)
    bc_bin = bc_axis.find(comp_frac)

    num_den_2d = numpy.zeros([diam_axis.n_bin, bc_axis.n_bin])
    diam_num_den = numpy.zeros([diam_axis.n_bin])
    diam_mass_den = numpy.zeros([diam_axis.n_bin])
    bc_num_den = numpy.zeros([bc_axis.n_bin])
    bc_mass_den = numpy.zeros([bc_axis.n_bin])

    for i in range(particles.n_particles):
        scale_2d = particles.comp_vol[i] * diam_axis.grid_size(diam_bin[i]) \
            * (bc_axis.grid_size(bc_bin[i]) / 100)
        diam_scale = particles.comp_vol[i] * diam_axis.grid_size(diam_bin[i])
        bc_scale = particles.comp_vol[i] * (bc_axis.grid_size(bc_bin[i]) / 100)
        
        if diam_axis.valid_bin(diam_bin[i]) and bc_axis.valid_bin(bc_bin[i]):
            num_den_2d[diam_bin[i]][bc_bin[i]] += 1.0 / scale_2d * 1e-6 # m^{-3} to cm^{-3}
        if diam_axis.valid_bin(diam_bin[i]):
            diam_num_den[diam_bin[i]] += 1.0 / diam_scale * 1e-6 # m^{-3} to cm^{-3}
            diam_mass_den[diam_bin[i]] += mass[i] / diam_scale
        if bc_axis.valid_bin(bc_bin[i]):
            bc_num_den[bc_bin[i]] += 1.0 / bc_scale * 1e-6 # m^{-3} to cm^{-3}
            bc_mass_den[bc_bin[i]] += mass[i] / bc_scale

    max_val = num_den_2d.max()
    value = num_den_2d / max_val
    plot_data_2d = pmc_histogram_2d_multi([value], diam_axis, bc_axis)
    diam_num_plot_data = [[diam_axis.center(i), diam_num_den[i]]
                          for i in range(diam_axis.n_bin)]
    diam_mass_plot_data = [[diam_axis.center(i), diam_mass_den[i]]
                           for i in range(diam_axis.n_bin)]
    bc_num_plot_data = [[bc_axis.center(i), bc_num_den[i]]
                        for i in range(bc_axis.n_bin)]
    bc_mass_plot_data = [[bc_axis.center(i), bc_mass_den[i]]
                         for i in range(bc_axis.n_bin)]

    return (plot_data_2d, max_val, diam_num_plot_data, diam_mass_plot_data,
            bc_num_plot_data, bc_mass_plot_data, env_state)

for use_coag in [True, False]:
    if use_coag:
        time_filename_list = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
    else:
        time_filename_list = get_time_filename_list(netcdf_dir_nc, netcdf_pattern_nc)
    for use_color in [True, False]:

        c = canvas.canvas()
        g21 = c.insert(graph.graphxy(
            width = grid_graph_width,
            x = graph.axis.log(min = diameter_axis_min,
                               max = diameter_axis_max,
                               title = diam_axis_label),
            y = graph.axis.linear(min = bc_axis_min,
                                  max = bc_axis_max,
                                  density = 1.2,
                                  title = bc_axis_label)))
        g11 = c.insert(graph.graphxy(
                width = grid_graph_width,
                ypos = g21.height + grid_h_space,
                x = graph.axis.linkedaxis(g21.axes["x"],
                                          painter = linked_grid_painter),
                y = graph.axis.log(min = diam_num_min_max[0],
                                   max = diam_num_min_max[1],
                                   title = diam_num_axis_label,
                                   painter = major_grid_painter),
                y2 = graph.axis.log(min = diam_mass_min_max[0],
                                    max = diam_mass_min_max[1],
                                    title = diam_mass_axis_label)))
        g22 = c.insert(graph.graphxy(
            width = grid_graph_width,
            xpos = g21.width + grid_h_space,
            x = graph.axis.log(min = bc_num_min_max[0],
                               max = bc_num_min_max[1],
                               title = bc_num_axis_label,
                               painter = major_grid_painter),
            x3 = graph.axis.log(min = bc_mass_min_max[0],
                                max = bc_mass_min_max[1],
                                title = bc_mass_axis_label),
            y = graph.axis.linkedaxis(g21.axes["y"],
                                      painter = linked_grid_painter)))

        g11.doaxes()
        g22.doaxes()

        time = time_hour * 3600.0
        filename = file_filename_at_time(time_filename_list, time)
        (plot_data_2d, max_val, diam_num_plot_data, diam_mass_plot_data,
         bc_num_plot_data, bc_mass_plot_data, env_state) = get_plot_data(filename)
        if use_color:
            attrs_num = [color_list[0], style.linewidth.Thick]
            attrs_mass = [color_list[1], style.linewidth.Thick]
            palette = rainbow_palette
        else:
            attrs_num = [line_style_list[0], style.linewidth.Thick]
            attrs_mass = [line_style_list[1], style.linewidth.Thick]
            palette = gray_palette

        g21.plot(graph.data.points(plot_data_2d,
                                   xmin = 1, xmax = 2, ymin = 3, ymax = 4,
                                   color = 5),
                 styles = [hsb_rect(palette)])
        g21.dolayout()
        for axisname in ["x", "y"]:
            for t in g21.axes[axisname].data.ticks:
                if t.ticklevel is not None:
                    g21.stroke(g21.axes[axisname].positioner.vgridpath(t.temp_v),
                               [style.linestyle.dotted])
        g21.dodata()
        g21.doaxes()

        g11.plot(graph.data.points(diam_num_plot_data,
                                 x = 1, y = 2),
               styles = [graph.style.line(lineattrs = attrs_num)])
        g11.plot(graph.data.points(diam_mass_plot_data,
                                 x = 1, y2 = 2),
               styles = [graph.style.line(lineattrs = attrs_mass)])

        g22.plot(graph.data.points(bc_num_plot_data,
                                 x = 2, y = 1),
               styles = [graph.style.line(lineattrs = attrs_num)])
        g22.plot(graph.data.points(bc_mass_plot_data,
                                 x3 = 2, y = 1),
               styles = [graph.style.line(lineattrs = attrs_mass)])

        if use_coag:
            extra_text = "with coagulation"
        else:
            extra_text = "no coagulation"
        write_time(g11, env_state)
        #boxed_text(g11, extra_text, point = [1, 1], anchor_point_rel = [1, 1])

        label_plot_line_boxed(g11, diam_num_plot_data, 0.08,
                        "number", [1, 1])
        label_plot_line_boxed(g11, diam_mass_plot_data, 0.05,
                        "mass", [0, 0], yaxis = g11.axes["y2"])
        label_plot_line_boxed(g22, bc_num_plot_data, 20,
                        "number", [1, 0], flip_xy = True)
        label_plot_line_boxed(g22, bc_mass_plot_data, 31,
                        "mass", [0, 1], flip_xy = True,
                        xaxis = g22.axes["x3"], label_offset = 0.7 * unit.v_mm)

        add_horiz_color_bar(g21,
                            min = 0.0,
                            max = max_val,
                            title = r"number conc. $n_{\rm BC,dry}(D,w)$ ($\rm cm^{-3}$)",
                            palette = palette,
                            bar_offset = 2.1,
                            above = False)

        if use_coag:
            coag_suffix = "wc"
        else:
            coag_suffix = "nc"
        if use_color:
            color_suffix = "color"
        else:
            color_suffix = "bw"
        out_filename = "%s_%s_%s.pdf" % (out_prefix, coag_suffix, color_suffix)
        c.writePDFfile(out_filename)
        if not use_color:
            print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
            print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
