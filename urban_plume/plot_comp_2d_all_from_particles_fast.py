#!/usr/bin/env python
# Copyright (C) 2007-2008 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, math
import copy as module_copy
from Scientific.IO.NetCDF import *
from pyx import *
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
import numpy

netcdf_dir = "out"
netcdf_re = re.compile(r"urban_plume_0.5_3am_wc_state_0001_([0-9]{8})\.nc")
netcdf_re = re.compile(r"urban_plume_0.5_3am_wc_state_0001_(00000201)\.nc")
netcdf_re = re.compile(r"urban_plume_0.5_3am_state_0001_([0-9]{7}1)\.nc")

show_particles = [[1029, [0.02, 70], "wet diesel"],
                  [1892, [0.02, 45], "wet gasoline"],
                  [4116, [0.5, 70], "dry diesel"],
                  #[4167, [0.5, 40], "dry gasoline"],
                  [4186, [0.8, 45], "dry gasoline"],
                  ]
show_particles = []

def process_file(in_filename, out_filename):
    ncf = NetCDFFile(in_filename)
    particles = aero_particle_array_t(ncf)
    x_axis = pmc_log_axis(min = 1e-2, max = 2, n_bin = 160)
    y_axis = pmc_linear_axis(min = 0, max = 100, n_bin = 100)
    water_frac_array = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
    mask_array = numpy.zeros([x_axis.n_bin, y_axis.n_bin], dtype=int)
    dry_num_den_array = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
    wet_num_den_array = numpy.zeros([x_axis.n_bin, y_axis.n_bin])
    diameter = particles.dry_diameter() * 1e6
    comp_frac = particles.mass(include = ["BC"]) \
                / particles.mass(exclude = ["H2O"]) * 100
    water_frac = particles.mass(include = ["H2O"]) \
                 / particles.mass() * 100
    #oc_frac = particles.mass(include = ["BC"]) \
    #          / particles.mass(include = ["BC", "OC"]) * 100
    x_bin = x_axis.find(diameter)
    y_bin = y_axis.find(comp_frac)
    show_coords = [[] for show in show_particles]
    for i in range(particles.n_particles):
        #bin_array[x_bin[i], y_bin[i]] += 1.0 / particles.comp_vol[i] \
        #                                 / x_axis.grid_size(x_bin[i]) \
        #                                 / y_axis.grid_size(y_bin[i])
        if water_frac[i] > 0:
            if water_frac_array[x_bin[i], y_bin[i]] > 0.0:
                water_frac_array[x_bin[i], y_bin[i]] \
                       = min(water_frac[i],
                             water_frac_array[x_bin[i], y_bin[i]])
            else:
                water_frac_array[x_bin[i], y_bin[i]] = water_frac[i]
            wet_num_den_array[x_bin[i], y_bin[i]] \
                                        += 1.0 / particles.comp_vol[i] \
                                        / x_axis.grid_size(x_bin[i]) \
                                        / y_axis.grid_size(y_bin[i])
        else:
            #dry_num_den_array[x_bin[i], y_bin[i]] += 1.0 / particles.comp_vol[i] \
            #                     / x_axis.grid_size(x_bin[i]) \
            #                     / y_axis.grid_size(y_bin[i])
            dry_num_den_array[x_bin[i], y_bin[i]] = 1.0
        for j, show in enumerate(show_particles):
            if particles.id[i] == show[0]:
                show_coords[j] = [diameter[i], comp_frac[i]]
        #if (oc_frac[i] > 10) and (water_frac[i] == 0.0):
        #    print "time =", ncf.variables["time"].getValue()
        #    print "i=%6d d=%10e b=%10e w=%10e o=%10e" \
        #          % (particles.id[i], diameter[i],
        #             comp_frac[i], water_frac[i], oc_frac[i])
    #max_val = bin_array.max()
    #max_val = 1e10
    max_val = 100
    water_frac_array = water_frac_array / max_val
    water_frac_array = water_frac_array.clip(0.0, 1.0)

    g = graph.graphxy(
        width = 8,
        x = graph.axis.log(min = x_axis.min,
                           max = x_axis.max,
                           title = r'dry diameter ($\mu$m)'),
        y = graph.axis.linear(min = y_axis.min,
                              max = y_axis.max,
                              title = r"$f_{\rm BC,all}$",
                              texter = graph.axis.texter.decimal(suffix
                                                                 = r"\%")))

    start_time = 3 * 3600.0
    time_of_day = start_time + float(ncf.variables["time"].getValue())
    time_of_day = time_of_day % (24 * 3600.0)
    hours = int(time_of_day / 3600.0)
    minutes = int(time_of_day / 60.0) % 60
    seconds = int(time_of_day) % 60

    if dry_num_den_array.max() > 0.0:
        g.plot(graph.data.points(pmc_histogram_2d(dry_num_den_array, x_axis, y_axis),
                                 xmin = 1, xmax = 2, ymin = 3, ymax = 4,
                                 color = 5),
               styles = [graph.style.rect(gray_palette)])


    saturation_min = 0.00001
    max_wet_num_den = wet_num_den_array.max()
    saturation = wet_num_den_array / max_wet_num_den
    #saturation = saturation ** (1.0/3.0)

    critical_point = 0.05
    sharpness = 2.0
    trans_power = log(0.5) / log(critical_point)
    trans_saturation = 2.0 * saturation**trans_power - 1.0
    sign_trans_saturation = where(trans_saturation >= 0.0,
                                  ones(shape(trans_saturation)),
                                  -ones(shape(trans_saturation)))
    saturation = ((abs(trans_saturation))**(1.0/sharpness)
                  * sign_trans_saturation + 1.0) / 2.0

    #saturation = saturation.clip(0.0, 1.0)
    #saturation = saturation * (1.0 - saturation_min) + saturation_min
    rects = pmc_histogram_2d_multi([water_frac_array, saturation],
                                    x_axis, y_axis)
    if water_frac_array.max() > 0.0:
        g.plot(graph.data.points(rects,
                                 xmin = 1, xmax = 2, ymin = 3, ymax = 4,
                                 color = 5, saturation = 6),
               styles = [hsb_rect(rainbow_palette)])
    boxed_text(g, 0.04, 0.9, "%02d:%02d LST" % (hours, minutes))
    for i in range(len(show_particles)):
        if len(show_coords[i]) > 0:
            label_point(g, show_coords[i][0], show_coords[i][1],
                        show_particles[i][1][0], show_particles[i][1][1],
                        show_particles[i][2])
    add_color_bar(g, min = 0.0, max = max_val,
                  title = r"water fraction", palette = rainbow_palette,
                  texter = graph.axis.texter.decimal(suffix = r"\%"),
                  bar_x_offset = 0.8)
    g.writePDFfile(out_filename)
    ncf.close()


filenames = os.listdir(netcdf_dir)
for filename in filenames:
    match = netcdf_re.search(filename)
    if match:
        output_key = match.group(1)
        in_filename = os.path.join(netcdf_dir, filename)
        out_filename = os.path.join(netcdf_dir,
                                    "comp_2d_bc_all_water_%s.pdf" % output_key)
        print in_filename
        process_file(in_filename, out_filename)
