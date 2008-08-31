#!/usr/bin/env python
# Copyright (C) 2007, 2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
from fig_helper import *

out_prefix = "figs/test_brownian"

brownian_netcdf_dir = r"../test/brownian/out"
brownian_netcdf_pattern = r"brown_mc_state_0001_([0-9]{8})\.nc"

times_hour = [0, 12, 24]

num_labels = [
    {"diameter": 0.08, "pos": [0, 1]},
    {"diameter": 0.17, "pos": [0, 1]},
    {"diameter": 0.27, "pos": [0, 1]},
    ]
              
vol_labels = [
    {"diameter": 0.08, "pos": [1, 1]},
    {"diameter": 0.17, "pos": [1, 1]},
    {"diameter": 0.22, "pos": [1, 1]},
    ]
              
x_axis = pmc_log_axis(min = diameter_axis_min, max = diameter_axis_max,
                      n_bin = num_diameter_bins)

times_sec = [t * 3600 for t in times_hour]

sect_data = pmc_var(NetCDFFile(os.path.join(brownian_netcdf_dir,
                                            "brown_sect_0001.nc")),
		    "aero",
		    [sum("aero_species")])

sect_data.scale_dim("radius",2e6)

sect_data.scale(math.log(10.0)) # d/dln(r) to d/dlog10(r)

for use_color in [True, False]:
    c = canvas.canvas()

    g_vol_lin = c.insert(graph.graphxy(
        width = 6.1,
        x = graph.axis.log(min = 1e-2,
                           max = 1e0,
                           title = r"dry diameter ($\rm \mu m$)",
                           painter = grid_painter),
        y = graph.axis.linear(min = 0,
                              max = 2.5e-7,
                              title = r"mass density ($\rm kg\,m^{-3}$)",
                              painter = grid_painter)))

    g_num_lin = c.insert(graph.graphxy(
        width = 6.1,
        ypos = g_vol_lin.height + 0.5,
        x = graph.axis.linkedaxis(g_vol_lin.axes["x"],
                                  painter = graph.axis.painter.linked(gridattrs = [style.linestyle.dotted])),
        y = graph.axis.linear(min = 0,
                              max = 2e11,
                              title = r"number density ($\rm m^{-3}$)",
                              painter = grid_painter)))

    time_filename_list = get_time_filename_list(brownian_netcdf_dir,
                                                brownian_netcdf_pattern)
    for i in range(len(times_sec)):
        filename = file_filename_at_time(time_filename_list, times_sec[i])
        ncf = NetCDFFile(filename)
        particles = aero_particle_array_t(ncf)
        ncf.close()

        diameter = particles.diameter() * 1e6
        mass = particles.mass() * 1.8 # scale density from H2O to 1800 kg/m^3

        x_bin = x_axis.find(diameter)

        num_den_array = numpy.zeros([x_axis.n_bin])
        mass_den_array = numpy.zeros([x_axis.n_bin])
        for k in range(particles.n_particles):
            scale = particles.comp_vol[k] * x_axis.grid_size(x_bin[k])
            num_den_array[x_bin[k]] += 1.0 / scale
            mass_den_array[x_bin[k]] += mass[k] / scale

        num_plot_data = [[x_axis.center(k), num_den_array[k]]
                         for k in range(x_axis.n_bin) if num_den_array[k] > 0.0]
        mass_plot_data = [[x_axis.center(k), mass_den_array[k]]
                         for k in range(x_axis.n_bin) if mass_den_array[k] > 0.0]

        if use_color:
            plot_color = color_list[i]
        else:
            plot_color = color.gray.black
        g_num_lin.plot(
            graph.data.points(num_plot_data, x = 1, y = 2,
                              title = "%g hours MC" % times_hour[i]),
            styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
                                         size = 0.08,
                                         symbolattrs = [plot_color])])

        g_vol_lin.plot(
            graph.data.points(mass_plot_data, x = 1, y = 2,
                              title = "%g hours MC" % times_hour[i]),
            styles = [graph.style.symbol(symbol = graph.style.symbol.circle,
                                         size = 0.08,
                                         symbolattrs = [plot_color])])

        data_slice = module_copy.deepcopy(sect_data)
        data_slice.reduce([select("unit", "num_den"),
                           select("time", times_sec[i])])

        plot_data = data_slice.data_center_list(strip_zero = True)
        g_num_lin.plot(
            graph.data.points(plot_data, x = 1, y = 2),
            styles = [graph.style.line(lineattrs = [plot_color])])
        label_plot_line(g_num_lin, plot_data, num_labels[i]["diameter"],
                        "%d hours" % times_hour[i],
                        label_pos = num_labels[i]["pos"])

        data_slice = module_copy.deepcopy(sect_data)
        data_slice.reduce([select("unit", "vol_den"),
                           select("time", times_sec[i])])
        data_slice.scale(1800) # volume to mass with density = 1800 kg/m^3
        plot_data = data_slice.data_center_list(strip_zero = True)
        g_vol_lin.plot(
            graph.data.points(plot_data, x = 1, y = 2),
            styles = [graph.style.line(lineattrs = [plot_color])])
        label_plot_line(g_vol_lin, plot_data, vol_labels[i]["diameter"],
                        "%d hours" % times_hour[i],
                        label_pos = vol_labels[i]["pos"])

    if use_color:
        out_filename = "%s_color.pdf" % out_prefix
    else:
        out_filename = "%s_bw.pdf" % out_prefix
    c.writePDFfile(out_filename)
    if use_color:
        print "figure height = %.1f cm" % unit.tocm(c.bbox().height())
        print "figure width = %.1f cm" % unit.tocm(c.bbox().width())
