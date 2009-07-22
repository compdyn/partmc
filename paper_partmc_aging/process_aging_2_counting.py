#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, pyx, re
sys.path.append("../tool")
from pmc_pyx import *
from numpy import *
sys.path.append("../tool")
from pmc_data_nc import *
from fig_helper import *

def delta(arr):
    return (arr[1:] - arr[:-1])

max_error_num = 0
max_error_mass = 0.0

for coag_suffix in ["wc", "nc"]:
    for (type_suffix, type) in [("num", int),
                                ("mass", float)]:
        print "Reading data %s %s..." % (coag_suffix, type_suffix)

        filename_pattern = "aging_%s_([0-9]{8})_time.txt" % coag_suffix
        filename_re = re.compile(filename_pattern)
        filename_list = get_filename_list(aging_data_dir, filename_pattern)
        n_time = len(filename_list)
        print "Found %d times" % n_time

        first_time = True
        for (i, filename) in enumerate(filename_list):
            print "Reading %s %s %d" % (coag_suffix, type_suffix, i)
            match = filename_re.search(filename)
            if not match:
                raise Exception()
            key = match.group(1)

            filename_base = os.path.join(aging_data_dir,
                                         "aging_%s_%s" % (coag_suffix, key))
            filename_base_type = os.path.join(aging_data_dir,
                                              "aging_%s_%s_%s" % (coag_suffix, key, type_suffix))

            time_array = loadtxt("%s_time.txt" % filename_base, float)
            height_array = loadtxt("%s_height.txt" % filename_base, float)
            comp_vol_array = loadtxt("%s_comp_vol.txt" % filename_base, float)

            data_total_array = loadtxt("%s.txt" % filename_base_type, type)
            data_emit_array = loadtxt("%s_emit.txt" % filename_base_type, type)
            data_dilution_array = loadtxt("%s_dilution.txt" % filename_base_type, type)
            data_halving_array = loadtxt("%s_halving.txt" % filename_base_type, type)
            data_cond_array = loadtxt("%s_cond.txt" % filename_base_type, type)
            data_coag_gain_array = loadtxt("%s_coag_gain.txt" % filename_base_type, type)
            data_coag_loss_array = loadtxt("%s_coag_loss.txt" % filename_base_type, type)

            if first_time:
                n_level = data_total_array.size - 1

                time = zeros([n_time], float)
                height = zeros([n_time], float)
                comp_vol = zeros([n_time], float)

                data_a = zeros([n_time, n_level], type)
                data_f = zeros([n_time, n_level], type)
                data_emit_a = zeros([n_time - 1, n_level], type)
                data_emit_f = zeros([n_time - 1, n_level], type)
                data_dilution_a = zeros([n_time - 1, n_level], type)
                data_dilution_f = zeros([n_time - 1, n_level], type)
                data_halving_a = zeros([n_time - 1, n_level], type)
                data_halving_f = zeros([n_time - 1, n_level], type)
                data_cond_a_a = zeros([n_time - 1, n_level], type)
                data_cond_a_f = zeros([n_time - 1, n_level], type)
                data_cond_f_a = zeros([n_time - 1, n_level], type)
                data_cond_f_f = zeros([n_time - 1, n_level], type)
                data_coag_gain_a = zeros([n_time - 1, n_level], type)
                data_coag_gain_f = zeros([n_time - 1, n_level], type)
                data_coag_loss_a_a = zeros([n_time - 1, n_level], type)
                data_coag_loss_a_f = zeros([n_time - 1, n_level], type)
                data_coag_loss_f_a = zeros([n_time - 1, n_level], type)
                data_coag_loss_f_f = zeros([n_time - 1, n_level], type)

            time[i] = float(time_array)
            height[i] = float(height_array)
            comp_vol[i] = float(comp_vol_array)

            for level in range(n_level):
                k = level + 1
                data_a[i, level] = data_total_array[:k].sum()
                data_f[i, level] = data_total_array[k:].sum()
                if not first_time:
                    data_emit_a[i - 1, level] = data_emit_array[:k].sum()
                    data_emit_f[i - 1, level] = data_emit_array[k:].sum()
                    data_dilution_a[i - 1, level] = data_dilution_array[:k].sum()
                    data_dilution_f[i - 1, level] = data_dilution_array[k:].sum()
                    data_halving_a[i - 1, level] = data_halving_array[:k].sum()
                    data_halving_f[i - 1, level] = data_halving_array[k:].sum()
                    data_cond_a_a[i - 1, level] = data_cond_array[:k,:k].sum()
                    data_cond_a_f[i - 1, level] = data_cond_array[:k,k:].sum()
                    data_cond_f_a[i - 1, level] = data_cond_array[k:,:k].sum()
                    data_cond_f_f[i - 1, level] = data_cond_array[k:,k:].sum()
                    data_coag_gain_a[i - 1, level] = data_coag_gain_array[:k].sum()
                    data_coag_gain_f[i - 1, level] = data_coag_gain_array[k:].sum()
                    data_coag_loss_a_a[i - 1, level] = data_coag_loss_array[:k,:k].sum()
                    data_coag_loss_a_f[i - 1, level] = data_coag_loss_array[:k,k:].sum()
                    data_coag_loss_f_a[i - 1, level] = data_coag_loss_array[k:,:k].sum()
                    data_coag_loss_f_f[i - 1, level] = data_coag_loss_array[k:,k:].sum()

            first_time = False

        for level in range(n_level):
            delta_data_a = delta(data_a[:,level])
            delta_data_f = delta(data_f[:,level])

            aged_kp1 = data_a[1:,level] - (data_cond_a_a[:,level] + data_cond_f_a[:,level] +
                                           data_emit_a[:,level] + data_coag_gain_a[:,level])
            aged_k = data_a[:-1,level] - (data_cond_a_a[:,level] + data_cond_a_f[:,level] +
                                   data_dilution_a[:,level] + data_halving_a[:,level] +
                                   data_coag_loss_a_a[:,level] + data_coag_loss_a_f[:,level])
            aged_delta = delta_data_a - (data_emit_a[:,level] - data_dilution_a[:,level] -
                                        data_halving_a[:,level] + data_cond_f_a[:,level] -
                                        data_cond_a_f[:,level] + data_coag_gain_a[:,level] -
                                        data_coag_loss_a_a[:,level] - data_coag_loss_a_f[:,level])
            fresh_kp1 = data_f[1:,level] - (data_cond_a_f[:,level] + data_cond_f_f[:,level] +
                                     data_emit_f[:,level] + data_coag_gain_f[:,level])
            fresh_k = data_f[:-1,level] - (data_cond_f_a[:,level] + data_cond_f_f[:,level] +
                                    data_dilution_f[:,level] + data_halving_f[:,level] +
                                    data_coag_loss_f_a[:,level] + data_coag_loss_f_f[:,level])
            fresh_delta = delta_data_f - (data_emit_f[:,level] - data_dilution_f[:,level] -
                                         data_halving_f[:,level] + data_cond_a_f[:,level] -
                                         data_cond_f_a[:,level] + data_coag_gain_f[:,level] -
                                         data_coag_loss_f_a[:,level] - data_coag_loss_f_f[:,level])

            max_error = max([max(abs(aged_kp1)), max(abs(aged_k)),
                             max(abs(aged_delta)), max(abs(fresh_kp1)),
                             max(abs(fresh_k)), max(abs(fresh_delta))])
            print "%s %4s %d -- %g" % (coag_suffix, type_suffix, level, max_error)
            if type_suffix == "num":
                max_error_num = max(max_error_num, max_error)
            else:
                max_error_mass = max(max_error_mass, max_error)

        filename = os.path.join(aging_data_dir,
                                "aging_%s_%s_%%s.txt" % (coag_suffix, type_suffix))
        if type_suffix == "num":
            fmt = "%d"
        else:
            fmt = "%.20e"
        savetxt(filename % "time", time, fmt = "%.20e")
        savetxt(filename % "height", height, fmt = "%.20e")
        savetxt(filename % "comp_vol", comp_vol, fmt = "%.20e")
        
        savetxt(filename % "a", data_a, fmt = fmt)
        savetxt(filename % "f", data_f, fmt = fmt)
        savetxt(filename % "emit_a", data_emit_a, fmt = fmt)
        savetxt(filename % "emit_f", data_emit_f, fmt = fmt)
        savetxt(filename % "dilution_a", data_dilution_a, fmt = fmt)
        savetxt(filename % "dilution_f", data_dilution_f, fmt = fmt)
        savetxt(filename % "halving_a", data_halving_a, fmt = fmt)
        savetxt(filename % "halving_f", data_halving_f, fmt = fmt)
        savetxt(filename % "cond_a_a", data_cond_a_a, fmt = fmt)
        savetxt(filename % "cond_a_f", data_cond_a_f, fmt = fmt)
        savetxt(filename % "cond_f_a", data_cond_f_a, fmt = fmt)
        savetxt(filename % "cond_f_f", data_cond_f_f, fmt = fmt)
        savetxt(filename % "coag_gain_a", data_coag_gain_a, fmt = fmt)
        savetxt(filename % "coag_gain_f", data_coag_gain_f, fmt = fmt)
        savetxt(filename % "coag_loss_a_a", data_coag_loss_a_a, fmt = fmt)
        savetxt(filename % "coag_loss_a_f", data_coag_loss_a_f, fmt = fmt)
        savetxt(filename % "coag_loss_f_a", data_coag_loss_f_a, fmt = fmt)
        savetxt(filename % "coag_loss_f_f", data_coag_loss_f_f, fmt = fmt)

print "max_error_num: ", max_error_num
print "max_error_mass: ", max_error_mass
