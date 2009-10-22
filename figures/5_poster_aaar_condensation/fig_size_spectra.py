#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math, re
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../../tool")
from pmc_data_nc import *
from pmc_pyx import *
import numpy
from config import *

out_prefix = "figs/size_spectra"

plot_times_secs = [0, 120, 600]

diameter_axis_max = 100.0

for [i_run, netcdf_pattern] in netcdf_indexed_patterns:
    out_filename = "%s_%d.pdf" % (out_prefix, i_run)
    print out_filename

    time_filename_list = get_time_filename_list(netcdf_dir, netcdf_pattern)
    filename = time_filename_list[0][1]
    ncf = NetCDFFile(filename)
    env_state = env_state_t(ncf)
    ncf.close()
    min_time_min = env_state.elapsed_time / 60.0

    g = graph.graphxy(
        width = graph_width,
        x = graph.axis.log(min = diameter_axis_min,
                           max = diameter_axis_max,
                           title = r'wet diameter $D\ (\rm\mu m)$',
                           painter = grid_painter),
        y = graph.axis.linear(min = 0,
                              max = 3e4,
                              title = r"number conc. ($\rm cm^{-3}$)",
                              painter = grid_painter),
        key = graph.key.key(pos = "tl",
                            keyattrs = [deco.stroked,
                                        deco.filled([color.rgb.white])],
                            hdist = 0.3 * unit.v_cm,
                            vdist = 0.2 * unit.v_cm))

    for (i_plot, plot_time_sec) in enumerate(plot_times_secs):
        time = plot_time_sec + min_time_min * 60.0
        filename = file_filename_at_time(time_filename_list, time)
        ncf = NetCDFFile(filename)
        particles = aero_particle_array_t(ncf)
        ncf.close()

        diameter = particles.diameter() * 1e6
        comp_frac = particles.mass(include = ["BC"]) \
                    / particles.mass(exclude = ["H2O"]) * 100

        x_axis = pmc_log_axis(min = diameter_axis_min, max = diameter_axis_max,
                              n_bin = num_diameter_bins)
        x_bin = x_axis.find(diameter)

        num_den_array = numpy.zeros([x_axis.n_bin])
        for i in range(particles.n_particles):
            if x_axis.valid_bin(x_bin[i]):
                scale = particles.comp_vol[i] * x_axis.grid_size(x_bin[i])
                num_den_array[x_bin[i]] += 1.0 / scale
        num_den_array /= 1e6

        plot_data = [[x_axis.center(i), num_den_array[i]]
                     for i in range(x_axis.n_bin)]
        attrs = [color_list[i_plot],
                 style.linewidth.normal]
        label = "%d min" % (plot_time_sec / 60)
        g.plot(graph.data.points(plot_data, x = 1, y = 2, title = label),
               styles = [graph.style.line(lineattrs = attrs)])

        if i_plot == len(plot_times_secs) - 1:
            cutoff = 2 # um
            tot_num_conc = 0.0
            active_num_conc = 0.0
            for i in range(len(diameter)):
                tot_num_conc += 1.0 / particles.comp_vol[i]
                if diameter[i] > cutoff:
                    active_num_conc += 1.0 / particles.comp_vol[i]
            tot_num_conc /= 1e6
            active_num_conc /= 1e6
            print "total number conc = %g cm^{-3}" % tot_num_conc
            print "activated number conc = %g cm^{-3}" % active_num_conc
            print "activated number frac = %g %%" % (active_num_conc / tot_num_conc * 100.0)

    #write_time(g, env_state)

    g.writePDFfile(out_filename)
    print "figure height = %.1f cm" % unit.tocm(g.bbox().height())
    print "figure width = %.1f cm" % unit.tocm(g.bbox().width())
