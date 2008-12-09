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

bc_fractions = [[0, 2],
                [2, 10],
                [60, 70],
                ]
small_diameter_limit = 0.05
large_diameter_limit = 0.1
base_small = []
base_large = []
new_small = []
new_large = []

out_prefix = "figs/aero_bc_mixing"

x_axis = pmc_log_axis(min = diameter_axis_min, max = diameter_axis_max,
                      n_bin = num_diameter_bins)

time_filename_list_wc = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
time_filename_list_nc = get_time_filename_list(netcdf_dir_nc, netcdf_pattern_nc)

for use_color in [True, False]:
    g = graph.graphxy(
            width = 6.7,
            x = graph.axis.log(min = x_axis.min,
                               max = x_axis.max,
                               title = r'dry diameter $D$ ($\rm \mu m$)',
                               painter = grid_painter),
            y = graph.axis.log(min = 1e-6,
                               max = 1e3,
                               title = r'mass conc. $m_{\rm BC,dry}(D,[w_1,w_2])$ ($\rm \mu g\, m^{-3}$)',
                               painter = grid_painter))

    for with_coag in [True, False]:
        if with_coag:
            filename = file_filename_at_time(time_filename_list_wc,
                                             time_hour * 3600)
        else:
            filename = file_filename_at_time(time_filename_list_nc,
                                             time_hour * 3600)
        ncf = NetCDFFile(filename)
        particles = aero_particle_array_t(ncf)
        env_state = env_state_t(ncf)
        ncf.close()

        diameter = particles.dry_diameter() * 1e6
        mass = particles.mass() * 1e9
        comp_frac = particles.mass(include = ["BC"]) \
                    / particles.mass(exclude = ["H2O"]) * 100
        x_bin = x_axis.find(diameter)

        for i_frac in range(len(bc_fractions)):
            mass_array = numpy.zeros([x_axis.n_bin])
            for i in range(particles.n_particles):
                if (comp_frac[i] >= bc_fractions[i_frac][0] * (1.0 - 1e-8)) \
                   and (comp_frac[i] <= bc_fractions[i_frac][1] * (1.0 + 1e8)):
                    scale = particles.comp_vol[i] * x_axis.grid_size(x_bin[i])
                    mass_array[x_bin[i]] += mass[i] / scale

            plot_data = [[x_axis.center(i), mass_array[i]]
                         for i in range(x_axis.n_bin)
                         if mass_array[i] > 0.0]

            attrs = []
            if with_coag:
                if use_color:
                    attrs.append(style.linestyle.solid)
                    attrs.append(style.linewidth.THick)
                else:
                    attrs.append(style.linewidth.THIck)
            else:
                if use_color:
                    attrs.append(style.linestyle.dashed)
                    attrs.append(style.linewidth.THick)
                else:
                    attrs.append(style.linewidth.Thick)
            if use_color:
                attrs.append(color_list[i_frac])
            else:
                attrs.append(line_style_list[i_frac])
            g.plot(graph.data.points(plot_data, x = 1, y = 2),
                   styles = [graph.style.line(lineattrs = attrs)])

            mass_den = mass / particles.comp_vol
            particle_select = (comp_frac >= bc_fractions[i_frac][0]) \
                              & (comp_frac <= bc_fractions[i_frac][1]) \
                              & (diameter < small_diameter_limit)
            val = (mass_den[particle_select]).sum()
            if not use_color:
                print ("coag = %s, BC frac %g%% to %g%%,"
                       " mass den below %g um = %g ug/m^3") \
                       % (str(with_coag), bc_fractions[i_frac][0],
                          bc_fractions[i_frac][1], small_diameter_limit, val)
                if with_coag:
                    new_small.append(val)
                else:
                    base_small.append(val)

            mass_den = mass / particles.comp_vol
            particle_select = (comp_frac >= bc_fractions[i_frac][0]) \
                              & (comp_frac <= bc_fractions[i_frac][1]) \
                              & (diameter > large_diameter_limit)
            val = (mass_den[particle_select]).sum()
            if not use_color:
                print ("coag = %s, BC frac %g%% to %g%%,"
                       " mass den above %g um = %g ug/m^3") \
                       % (str(with_coag), bc_fractions[i_frac][0],
                          bc_fractions[i_frac][1], large_diameter_limit, val)
                if with_coag:
                    new_large.append(val)
                else:
                    base_large.append(val)

    if not use_color:
        if (len(base_small) != len(bc_fractions)) \
                or (len(base_large) != len(bc_fractions)) \
                or (len(new_small) != len(bc_fractions)) \
                or (len(new_large) != len(bc_fractions)):
            print "ERROR: problem computing base/new_large/small"
        for i in range(len(bc_fractions)):
            print ("BC frac %g%% to %g%%, decrease in mass den below %g um"
                   " from no-coag to with-coag = %g%%") \
                   % (bc_fractions[i][0], bc_fractions[i][1],
                      small_diameter_limit,
                      (base_small[i] - new_small[i]) / base_small[i] * 100)
            print ("BC frac %g%% to %g%%, decrease in mass den above %g um"
                   " from no-coag to with-coag = %g%%") \
                   % (bc_fractions[i][0], bc_fractions[i][1],
                      large_diameter_limit,
                      (base_large[i] - new_large[i]) / base_large[i] * 100)

    g.doaxes()
    g.dodata()

    ######################################################################
    # key

    dx = 1 * unit.v_cm
    dy = 0.5 * unit.v_cm
    xoff = 0.7 * unit.v_cm
    yoff = 0.3 * unit.v_cm
    length = 0.6 * unit.v_cm
    extra = 0.2 * unit.v_cm
    coag_off = 0.4 * unit.v_cm
    cornerx = 0.3 * unit.v_cm
    cornery = 0.2 * unit.v_cm

    c = canvas.canvas()

    for i in range(len(bc_fractions)):
        p = c.text(0, - i * dy,
                   "%d--%d\\%% BC" % (bc_fractions[i][0], bc_fractions[i][1]),
                   [text.halign.right, text.valign.middle])
        left = p.bbox().left()
        bottom = p.bbox().bottom()
        if i == 0:
            first_top = p.bbox().top()

    with_text = c.text(xoff, yoff,
                       "with", [text.halign.center])
    without_text = c.text(xoff + dx, yoff,
                          "without", [text.halign.center])
    centerx = (with_text.bbox().left() + without_text.bbox().right()) / 2.0
    coag_text = c.text(centerx, yoff + coag_off,
                       "coagulation", [text.halign.center])

    for i in range(len(bc_fractions)):
        for with_coag in [True, False]:
            if with_coag:
                if use_color:
                    attrs.append(style.linestyle.solid)
                    attrs.append(style.linewidth.THick)
                else:
                    attrs.append(style.linewidth.THick)
            else:
                if use_color:
                    attrs.append(style.linestyle.dashed)
                    attrs.append(style.linewidth.THick)
                else:
                    attrs.append(style.linewidth.thick)
            if use_color:
                attrs.append(color_list[i])
            else:
                attrs.append(line_style_list[i])
            if with_coag:
                c.stroke(path.line(xoff - length / 2.0, - i * dy,
                                   xoff + length / 2.0, - i * dy),
                         attrs)
            else:
                c.stroke(path.line(xoff + dx - length / 2.0, - i * dy,
                                   xoff + dx + length / 2.0, - i * dy),
                         attrs)

    #box = path.rect(left - extra, bottom - extra,
    #                without_text.bbox().right() - left + 2.0 * extra,
    #                coag_text.bbox().top() - bottom + 2.0 * extra)
    box = path.path(path.moveto(left - extra, bottom - extra),
                    path.lineto(left - extra, first_top + extra),
                    path.lineto(with_text.bbox().left() - extra,
                                first_top + extra),
                    path.lineto(with_text.bbox().left() - extra,
                                coag_text.bbox().top() + extra),
                    path.lineto(without_text.bbox().right() + extra,
                                coag_text.bbox().top() + extra),
                    path.lineto(without_text.bbox().right() + extra,
                                bottom - extra),
                    path.closepath())
    desiredx = g.width - cornerx
    desiredy = cornery
    transx = - box.bbox().right() + desiredx
    transy = - box.bbox().bottom() + desiredy
    g.draw(box, [deco.stroked, deco.filled([color.gray.white]),
                 trafo.translate(transx, transy)])
    g.insert(c, [trafo.translate(transx, transy)])

    ######################################################################

    write_time_inside(g, env_state)

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    g.writePDFfile(out_filename)
    if not use_color:
        print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
