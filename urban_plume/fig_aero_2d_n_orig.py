#!/usr/bin/env python
# Copyright (C) 2007-2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
sys.path.append(".")
from fig_helper import *

time_hour = 24
max_n_coags = 20

out_prefix = "figs/aero_2d_n_orig"

time_filename_list = get_time_filename_list(netcdf_dir_wc, netcdf_pattern_wc)
for color in [True, False]:
    g = graph.graphxy(
        width = 7.1,
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
                           title = diameter_axis_label),
        y = graph.axis.linear(min = -0.5,
                              max = max_n_coags + 0.5,
                              parter = graph.axis.parter.linear(tickdists = [4, 2]),
                              title = 'number of coagulation events $k$'))

    time = time_hour * 3600.0
    filename = file_filename_at_time(time_filename_list, time)
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    env_state = env_state_t(ncf)
    ncf.close()

    diameter = particles.dry_diameter() * 1e6

    x_axis = pmc_log_axis(min = diameter_axis_min, max = diameter_axis_max,
                          n_bin = num_diameter_bins)
    y_axis = pmc_linear_axis(min = -0.5, max = max_n_coags + 0.5,
                             n_bin = max_n_coags + 1)
    x_bin = x_axis.find(diameter)

    num_den_array = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
    for i in range(particles.n_particles):
        scale = particles.comp_vol[i] * x_axis.grid_size(x_bin[i])
        n_coags = int(particles.n_orig_part[i]) - 1
        n_coags = min(n_coags, max_n_coags)
        num_den_array[x_bin[i], n_coags] += 1.0 / scale * 1e-6 # m^{-3} to cm^{-3}

    n_coags = array(particles.n_orig_part) - 1
    lower_n_coag = 5
    if not color:
        print "%g%% particles have at least %d coagulation events" \
              % (len(n_coags[n_coags >= lower_n_coag])
                 / float(len(n_coags)) * 100,
                 lower_n_coag)
        print "max coagulation events = %d" % n_coags.max()
        for d in [0.1, 0.3]:
            x_bin = x_axis.find(array([d]))[0]
            n_coag_list = [c for c in range(max_n_coags)
                           if num_den_array[x_bin, c] > 0]
            print "diameter %g um, num coags from %d to %d" \
                  % (d, min(n_coag_list), max(n_coag_list))

    value = num_den_array
    value_max = value.max()
    if value_max > 0.0:
        value = value / value_max
    value = value.clip(0.0, 1.0)

    rects = pmc_histogram_2d_multi([value],
                                   x_axis, y_axis)
    if color:
        palette = rainbow_palette
    else:
        palette = gray_palette
    g.plot(graph.data.points(rects,
                             xmin = 1, xmax = 2, ymin = 3, ymax = 4,
                             color = 5),
           styles = [hsb_rect(palette)])

    g.dolayout()
    for axisname in ["x", "y"]:
        for t in g.axes[axisname].data.ticks:
            if t.ticklevel is not None:
                g.stroke(g.axes[axisname].positioner.vgridpath(t.temp_v),
                         [style.linestyle.dotted])
    g.dodata()
    g.doaxes()

    write_time_inside(g, env_state)

    add_horiz_color_bar(g,
                        min = 0.0,
                        max = value_max,
                        title = r"number conc. $n_{\rm coag}(D,k)\ (\rm cm^{-3})$",
                        palette = palette,
                        bar_offset = 0.6)

    if color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    g.writePDFfile(out_filename)
    if not color:
        print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
