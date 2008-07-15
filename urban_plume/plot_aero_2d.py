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
netcdf_re = re.compile(r"urban_plume_state_0001_([0-9]{7}1)\.nc")
#netcdf_re = re.compile(r"urban_plume_state_0001_(00000901)\.nc")

show_particles = [[1049, [0.02, 75], "wet diesel"],
                  [1279, [0.02, 15], "wet gasoline"],
                  [2329, [0.7, 75], "dry diesel"],
                  [2315, [0.7, 15], "dry gasoline"],
                  ]
show_particles = []

def ratio_of_sum(array_1, array_2):
    array_sum = array_1 + array_2
    ratio = where(array_sum != 0.0,
                  array_1 / array_sum,
                  zeros(shape(array_sum)))
    return ratio

def sharpen_filter(x, sharpness, critical_point):
    """x must be an array with values between 0 and 1
    sharpness of 1 means no sharpening
    sharpness > 1 means increasingly sharper
    sharpness < 1 means less sharp
    critical point in [0,1] is inflection point of sharpener
    the output y = f(x) will satisfy:
    f(0) = 0
    f(critical_point) = 0.5
    f(1) = 1"""
    trans_power = log(0.5) / log(critical_point)
    trans_x = 2.0 * x**trans_power - 1.0
    sign_trans_x = where(trans_x >= 0.0,
                         ones(shape(trans_x)),
                         -ones(shape(trans_x)))
    y = ((abs(trans_x))**(1.0/sharpness) * sign_trans_x + 1.0) / 2.0
    return y

def white_black_dist_to_sb(dist_from_white, color_black_ratio):
    # (saturation, brightness) =
    #            (0, 1 - dist_from_white)  when color_black_ratio = 0
    #            (dist_from_white, 1)      when color_black_ratio = 1
    # interpolate linearly with color_black_ratio
    saturation = color_black_ratio * dist_from_white
    brightness = (1 - color_black_ratio) * (1 - dist_from_white) \
                 + color_black_ratio 
    return (saturation, brightness)

def process_file(in_filename, out_filename):
    ncf = NetCDFFile(in_filename)
    particles = aero_particle_array_t(ncf)
    x_axis = pmc_log_axis(min = 1e-2, max = 2, n_bin = 70)
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
        mask_array[x_bin[i],y_bin[i]] = 1
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
            dry_num_den_array[x_bin[i], y_bin[i]] \
                                        += 1.0 / particles.comp_vol[i] \
                                        / x_axis.grid_size(x_bin[i]) \
                                        / y_axis.grid_size(y_bin[i])
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
    #max_val = 100
    #water_frac_array = water_frac_array / max_val
    #water_frac_array = water_frac_array.clip(0.0, 1.0)

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

    start_time = float(ncf.variables["start_time_of_day"].getValue())
    time_of_day = start_time + float(ncf.variables["elapsed_time"].getValue())
    time_of_day = time_of_day % (24 * 3600.0)
    hours = int(time_of_day / 3600.0)
    minutes = int(time_of_day / 60.0) % 60
    seconds = int(time_of_day) % 60

    #if dry_num_den_array.max() > 0.0:
    #    g.plot(graph.data.points(pmc_histogram_2d(dry_num_den_array, x_axis, y_axis),
    #                             xmin = 1, xmax = 2, ymin = 3, ymax = 4,
    #                             color = 5),
    #           styles = [graph.style.rect(gray_palette)])


    #saturation_min = 0.00001
    #max_wet_num_den = wet_num_den_array.max()
    #saturation = wet_num_den_array / max_wet_num_den
    #saturation = saturation ** (1.0/3.0)

    #saturation = sharpen_filter(saturation,
    #                            sharpness = 0.2,
    #                            critical_point = 0.05)

    dist_from_white = dry_num_den_array + wet_num_den_array
    dist_from_white = dist_from_white / dist_from_white.max()
    dist_from_white.clip(0.0, 1.0)
    dist_from_white = sharpen_filter(dist_from_white,
                                sharpness = 0.2,
                                critical_point = 0.05)

    color_black_ratio = ratio_of_sum(wet_num_den_array, dry_num_den_array)

    (saturation, brightness) = white_black_dist_to_sb(dist_from_white,
                                                      color_black_ratio)

    hue = water_frac_array
    hue_max = hue.max()
    if hue_max > 0.0:
        hue = hue / hue.max()
    hue.clip(0.0, 1.0)

    hue = dry_num_den_array + wet_num_den_array
    hue = log(1 + hue / hue.max() * 100.0)
    #hue = hue / hue.sum() / x_axis.grid_size(0) \
    #        / (y_axis.grid_size(0) / 100.0)
    #hue_max = 2
    hue_max = hue.max()
    if hue_max > 0.0:
        hue = hue / hue_max
    hue = hue.clip(0.0, 1.0)

    saturation = ones(shape(hue))
    brightness = ones(shape(hue))
    #brightness = where(water_frac_array > 0.0,
    #              ones(shape(hue)),
    #              zeros(shape(hue)))

    #print "dry_num_den_array"
    #print dry_num_den_array[39,80:]
    #print "wet_num_den_array"
    #print wet_num_den_array[39,80:]
    #print "dist_from_white"
    #print dist_from_white[39,80:]
    #print "color_black_ratio"
    #print color_black_ratio[39,80:]
    #print "hue"
    #print hue[39,80:]
    #print "saturation"
    #print saturation[39,80:]
    #print "brightness"
    #print brightness[39,80:]
    
    #saturation = ones(shape(hue))
    #brightness = ones(shape(hue))

    #saturation = saturation.clip(0.0, 1.0)
    #saturation = saturation * (1.0 - saturation_min) + saturation_min
    rects = pmc_histogram_2d_multi([hue, saturation, brightness],
                                    x_axis, y_axis, mask = mask_array)
    #if water_frac_array.max() > 0.0:
    g.plot(graph.data.points(rects,
                             xmin = 1, xmax = 2, ymin = 3, ymax = 4,
                             color = 5, saturation = 6, brightness = 7),
           styles = [hsb_rect(rainbow_palette)])
    g.dodata()
    g.doaxes()

    boxed_text(g, 0.04, 0.9, "%02d:%02d LST" % (hours, minutes))
    for i in range(len(show_particles)):
        if len(show_coords[i]) > 0:
            label_point(g, show_coords[i][0], show_coords[i][1],
                        show_particles[i][1][0], show_particles[i][1][1],
                        show_particles[i][2])
    add_color_bar(g, min = 0.0, max = 100,
                  title = r"water fraction", palette = rainbow_palette,
                  texter = graph.axis.texter.decimal(suffix = r"\%"),
                  bar_x_offset = 0.8)

    g.fill(path.rect(g.width + 0.8, -0.5, 0.5, 0.5),
           [color.rgb.black])
    g.text(g.width + 0.8 + 0.5 + 0.3, -0.25, "dry",
           [text.halign.left, text.valign.middle])
    
    g.writePDFfile(out_filename)
    ncf.close()

def find_single_particles(early_filename, late_filename):
    early_ncf = NetCDFFile(early_filename)
    late_ncf = NetCDFFile(late_filename)
    early_particles = aero_particle_array_t(early_ncf)
    late_particles = aero_particle_array_t(late_ncf)
    early_diameter = early_particles.dry_diameter() * 1e6
    early_comp_frac = early_particles.mass(include = ["BC"]) \
                      / early_particles.mass(exclude = ["H2O"]) * 100
    early_water_frac = early_particles.mass(include = ["H2O"]) \
                       / early_particles.mass() * 100
    early_oc_frac = early_particles.mass(include = ["BC"]) \
                    / early_particles.mass(include = ["BC", "OC"]) * 100
    found = {}
    for i in range(early_particles.n_particles):
        if early_particles.id[i] in late_particles.id:
            type = ""
            if early_comp_frac[i] > 50:
                if early_water_frac[i] == 0.0:
                    type = "dry diesel"
                else:
                    type = "wet diesel"
            elif early_comp_frac[i] > 10:
                if early_water_frac[i] == 0.0:
                    type = "dry gasoline"
                else:
                    type = "wet gasoline"
            if type != "":
                if type not in found.keys():
                    found[type] = []
                found[type].append(i)
    for type in found.keys():
        found[type].sort(cmp = lambda x,y : cmp(early_particles.id[x],
                                                early_particles.id[y]))
        for i in found[type]:
            print "%s: id = %d, d = %f, bc = %g, water = %f" \
                  % (type, early_particles.id[i], early_diameter[i],
                     early_comp_frac[i], early_water_frac[i])

def find_closest():
    filename = "out/urban_plume_state_0001_00000781.nc"
    print "filename = %s" % filename
    ncf = NetCDFFile(filename)
    particles = aero_particle_array_t(ncf)
    diameter = particles.dry_diameter() * 1e6
    comp_frac = particles.mass(include = ["BC"]) \
                / particles.mass(exclude = ["H2O"]) * 100

    env_state = env_state_t(ncf)
    print "elapsed time = %g seconds" % env_state.elapsed_time
    print "time of day = %s LST" \
          % time_of_day_string(env_state.elapsed_time
                               + env_state.start_time_of_day)
    print "RH = %g %%" % (env_state.relative_humidity * 100)
    print "temp = %g K" % env_state.temperature
    print "pressure = %g Pa" % env_state.pressure
    print "height = %g m" % env_state.height

    gas_state = gas_state_t(ncf)
    print "gas concentrations:"
    for i in range(len(gas_state.concentration)):
        print "   %s = %g ppb" % (gas_state.gas_data.name[i],
                                  gas_state.concentration[i])

    for desired_diameter in [0.05, 0.07, 0.1]:
        for desired_comp_frac in [25, 50, 75]:
            i = ((diameter / desired_diameter - 1)**2
                 + (comp_frac / desired_comp_frac - 1)**2).argmin()
            print "**********************************************************"
            print ("searching for particle with diameter near %g um and "
                   "comp_frac near %g %%") \
                   % (desired_diameter, desired_comp_frac)
            print "found particle ID = %d" % particles.id[i]
            print "diameter = %g um" % diameter[i]
            print "comp_frac = %g %%" % comp_frac[i]
            print "consituent masses:"
            for j in range(len(particles.aero_data.name)):
                print "   %s = %g kg" % (particles.aero_data.name[j],
                                         particles.masses[j,i])

def make_figs():
    filenames = os.listdir(netcdf_dir)
    for filename in filenames:
        match = netcdf_re.search(filename)
        if match:
            output_key = match.group(1)
            in_filename = os.path.join(netcdf_dir, filename)
            out_filename = os.path.join(netcdf_dir,
                                        "comp_2d_bc_all_num_%s.pdf"
                                        % output_key)
            print in_filename
            process_file(in_filename, out_filename)

make_figs()

#find_single_particles("out/urban_plume_state_0001_00000301.nc", # early
#                      "out/urban_plume_state_0001_00001441.nc") # late

#find_closest()
