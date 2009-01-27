#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, pyx
sys.path.append("../tool")
from pmc_data_nc import *
from pmc_pyx import *
from numpy import *

def mean_day(time, data):
    return mean([x for [t,x] in filter_inf(zip(time, data))
                 if t >= 6 * 3600 - 0.1 and t <= 9 * 3600 + 0.1])

def mean_night(time, data):
    return mean([x for [t,x] in filter_inf(zip(time, data))
                 if t >= 12 * 3600 - 0.1 and t <= 22 * 3600 + 0.1])

dilution_rate = 1.5e-5 # s^{-1}

data_prefix = "out"

smooth_window_len = 60

grey_level = 0.2

max_error_num = 0
max_error_mass = 0.0

for coag_suffix in ["wc", "nc"]:
    for (type_suffix, type) in [("num", int),
                                ("mass", float)]:
        print "Reading data %s %s..." % (coag_suffix, type_suffix)

        filename = "%s/aging_%s_%s_%%s.txt" % (data_prefix, coag_suffix, type_suffix)
        if type_suffix == "num":
            type = int
        else:
            type = float

        time = loadtxt(filename % "time", float)
        height = loadtxt(filename % "height", float)
        comp_vol = loadtxt(filename % "comp_vol", float)

        data_a = loadtxt(filename % "a", type)
        data_f = loadtxt(filename % "f", type)
        data_emit_a = loadtxt(filename % "emit_a", type)
        data_emit_f = loadtxt(filename % "emit_f", type)
        data_dilution_a = loadtxt(filename % "dilution_a", type)
        data_dilution_f = loadtxt(filename % "dilution_f", type)
        data_halving_a = loadtxt(filename % "halving_a", type)
        data_halving_f = loadtxt(filename % "halving_f", type)
        data_cond_a_a = loadtxt(filename % "cond_a_a", type)
        data_cond_a_f = loadtxt(filename % "cond_a_f", type)
        data_cond_f_a = loadtxt(filename % "cond_f_a", type)
        data_cond_f_f = loadtxt(filename % "cond_f_f", type)
        data_coag_gain_a = loadtxt(filename % "coag_gain_a", type)
        data_coag_gain_f = loadtxt(filename % "coag_gain_f", type)
        data_coag_loss_a_a = loadtxt(filename % "coag_loss_a_a", type)
        data_coag_loss_a_f = loadtxt(filename % "coag_loss_a_f", type)
        data_coag_loss_f_a = loadtxt(filename % "coag_loss_f_a", type)
        data_coag_loss_f_f = loadtxt(filename % "coag_loss_f_f", type)

        n_time = time.shape[0]
        n_level = data_a.shape[1]

        max_time_min = max(time) / 60
        start_time_of_day_min = 6 * 60

        data_a_smooth = zeros([n_time, n_level], type)
        data_f_smooth = zeros([n_time, n_level], type)
        data_emit_a_smooth = zeros([n_time - 1, n_level], type)
        data_emit_f_smooth = zeros([n_time - 1, n_level], type)
        data_dilution_a_smooth = zeros([n_time - 1, n_level], type)
        data_dilution_f_smooth = zeros([n_time - 1, n_level], type)
        data_halving_a_smooth = zeros([n_time - 1, n_level], type)
        data_halving_f_smooth = zeros([n_time - 1, n_level], type)
        data_cond_a_a_smooth = zeros([n_time - 1, n_level], type)
        data_cond_a_f_smooth = zeros([n_time - 1, n_level], type)
        data_cond_f_a_smooth = zeros([n_time - 1, n_level], type)
        data_cond_f_f_smooth = zeros([n_time - 1, n_level], type)
        data_coag_gain_a_smooth = zeros([n_time - 1, n_level], type)
        data_coag_gain_f_smooth = zeros([n_time - 1, n_level], type)
        data_coag_loss_a_a_smooth = zeros([n_time - 1, n_level], type)
        data_coag_loss_a_f_smooth = zeros([n_time - 1, n_level], type)
        data_coag_loss_f_a_smooth = zeros([n_time - 1, n_level], type)
        data_coag_loss_f_f_smooth = zeros([n_time - 1, n_level], type)

        k_transfer_cond = zeros([n_time - 1, n_level], float)
        k_transfer_cond_net = zeros([n_time - 1, n_level], float)
        k_transfer = zeros([n_time - 1, n_level], float)
        k_transfer_net = zeros([n_time - 1, n_level], float)

        k_transfer_cond_smooth = zeros([n_time - 1, n_level], float)
        k_transfer_cond_net_smooth = zeros([n_time - 1, n_level], float)
        k_transfer_smooth = zeros([n_time - 1, n_level], float)
        k_transfer_net_smooth = zeros([n_time - 1, n_level], float)

        tau_transfer_cond = zeros([n_time - 1, n_level], float)
        tau_transfer_cond_net = zeros([n_time - 1, n_level], float)
        tau_transfer = zeros([n_time - 1, n_level], float)
        tau_transfer_net = zeros([n_time - 1, n_level], float)

        tau_transfer_cond_smooth = zeros([n_time - 1, n_level], float)
        tau_transfer_cond_net_smooth = zeros([n_time - 1, n_level], float)
        tau_transfer_smooth = zeros([n_time - 1, n_level], float)
        tau_transfer_net_smooth = zeros([n_time - 1, n_level], float)

        tau_day = zeros([n_level], float)
        tau_day_cond = zeros([n_level], float)
        tau_night = zeros([n_level], float)
        tau_night_cond = zeros([n_level], float)

        for level in range(n_level):
            data_a_smooth[:,level] = smooth(data_a[:,level], window_len = smooth_window_len)
            data_f_smooth[:,level] = smooth(data_f[:,level], window_len = smooth_window_len)
            data_emit_a_smooth[:,level] = smooth(data_emit_a[:,level], window_len = smooth_window_len)
            data_emit_f_smooth[:,level] = smooth(data_emit_f[:,level], window_len = smooth_window_len)
            data_dilution_a_smooth[:,level] = smooth(data_dilution_a[:,level], window_len = smooth_window_len)
            data_dilution_f_smooth[:,level] = smooth(data_dilution_f[:,level], window_len = smooth_window_len)
            data_halving_a_smooth[:,level] = smooth(data_halving_a[:,level], window_len = smooth_window_len)
            data_halving_f_smooth[:,level] = smooth(data_halving_f[:,level], window_len = smooth_window_len)
            data_cond_a_a_smooth[:,level] = smooth(data_cond_a_a[:,level], window_len = smooth_window_len)
            data_cond_a_f_smooth[:,level] = smooth(data_cond_a_f[:,level], window_len = smooth_window_len)
            data_cond_f_a_smooth[:,level] = smooth(data_cond_f_a[:,level], window_len = smooth_window_len)
            data_cond_f_f_smooth[:,level] = smooth(data_cond_f_f[:,level], window_len = smooth_window_len)
            data_coag_gain_a_smooth[:,level] = smooth(data_coag_gain_a[:,level], window_len = smooth_window_len)
            data_coag_gain_f_smooth[:,level] = smooth(data_coag_gain_f[:,level], window_len = smooth_window_len)
            data_coag_loss_a_a_smooth[:,level] = smooth(data_coag_loss_a_a[:,level], window_len = smooth_window_len)
            data_coag_loss_a_f_smooth[:,level] = smooth(data_coag_loss_a_f[:,level], window_len = smooth_window_len)
            data_coag_loss_f_a_smooth[:,level] = smooth(data_coag_loss_f_a[:,level], window_len = smooth_window_len)
            data_coag_loss_f_f_smooth[:,level] = smooth(data_coag_loss_f_f[:,level], window_len = smooth_window_len)

            data_delta_a = delta(data_a[:,level])
            data_delta_f = delta(data_f[:,level])

            aged_kp1 = data_a[1:,level] - (data_cond_a_a[:,level] + data_cond_f_a[:,level] +
                                    data_emit_a[:,level] + data_coag_gain_a[:,level])
            aged_k = data_a[:-1,level] - (data_cond_a_a[:,level] + data_cond_a_f[:,level] +
                                   data_dilution_a[:,level] + data_halving_a[:,level] +
                                   data_coag_loss_a_a[:,level] + data_coag_loss_a_f[:,level])
            aged_delta = data_delta_a - (data_emit_a[:,level] - data_dilution_a[:,level] -
                                        data_halving_a[:,level] + data_cond_f_a[:,level] -
                                        data_cond_a_f[:,level] + data_coag_gain_a[:,level] -
                                        data_coag_loss_a_a[:,level] - data_coag_loss_a_f[:,level])
            fresh_kp1 = data_f[1:,level] - (data_cond_a_f[:,level] + data_cond_f_f[:,level] +
                                     data_emit_f[:,level] + data_coag_gain_f[:,level])
            fresh_k = data_f[:-1,level] - (data_cond_f_a[:,level] + data_cond_f_f[:,level] +
                                    data_dilution_f[:,level] + data_halving_f[:,level] +
                                    data_coag_loss_f_a[:,level] + data_coag_loss_f_f[:,level])
            fresh_delta = data_delta_f - (data_emit_f[:,level] - data_dilution_f[:,level] -
                                         data_halving_f[:,level] + data_cond_a_f[:,level] -
                                         data_cond_f_a[:,level] + data_coag_gain_f[:,level] -
                                         data_coag_loss_f_a[:,level] - data_coag_loss_f_f[:,level])

            #print "aged k+1: ", aged_kp1
            #print "aged k: ", aged_k
            #print "aged delta: ", aged_delta
            #print "fresh k+1: ", fresh_kp1
            #print "fresh k: ", fresh_k
            #print "fresh delta: ", fresh_delta

            #print "max(aged k+1): ", max(abs(aged_kp1))
            #print "max(aged k): ", max(abs(aged_k))
            #print "max(aged delta): ", max(abs(aged_delta))
            #print "max(fresh k+1): ", max(abs(fresh_kp1))
            #print "max(fresh k): ", max(abs(fresh_k))
            #print "max(fresh delta): ", max(abs(fresh_delta))

            max_error = max([max(abs(aged_kp1)), max(abs(aged_k)),
                max(abs(aged_delta)), max(abs(fresh_kp1)),
                max(abs(fresh_k)), max(abs(fresh_delta))])
            print "%s %4s %d -- %g" % (coag_suffix, type_suffix, level, max_error)
            if type_suffix == "num":
                max_error_num = max(max_error_num, max_error)
            else:
                max_error_mass = max(max_error_mass, max_error)

            k_transfer_cond[:,level] = data_cond_f_a[:,level] / delta(time) / data_f[1:,level]
            k_transfer_cond_net[:,level] = (data_cond_f_a[:,level] - data_cond_a_f[:,level]) / delta(time) / data_f[1:,level]
            k_transfer[:,level] = (data_cond_f_a[:,level] + data_coag_loss_f_a[:,level]) / delta(time) / data_f[1:,level]
            k_transfer_net[:,level] = (data_cond_f_a[:,level] - data_cond_a_f[:,level]
                              + data_coag_loss_f_a[:,level] - data_coag_loss_a_f[:,level]) \
                              / delta(time) / data_f[1:,level]

            k_transfer_cond_smooth[:,level] = data_cond_f_a_smooth[:,level] / delta(time) / data_f_smooth[1:,level]
            k_transfer_cond_net_smooth[:,level] = (data_cond_f_a_smooth[:,level] - data_cond_a_f_smooth[:,level]) / delta(time) / data_f_smooth[1:,level]
            k_transfer_smooth[:,level] = (data_cond_f_a_smooth[:,level] + data_coag_loss_f_a_smooth[:,level]) / delta(time) / data_f_smooth[1:,level]
            k_transfer_net_smooth[:,level] = (data_cond_f_a_smooth[:,level] - data_cond_a_f_smooth[:,level]
                                     + data_coag_loss_f_a_smooth[:,level] - data_coag_loss_a_f_smooth[:,level]) \
                                     / delta(time) / data_f_smooth[1:,level]

            k_transfer_cond[:,level] = nan_to_value(k_transfer_cond[:,level], 0.0)
            k_transfer_cond_net[:,level] = nan_to_value(k_transfer_cond_net[:,level], 0.0)
            k_transfer[:,level] = nan_to_value(k_transfer[:,level], 0.0)
            k_transfer_net[:,level] = nan_to_value(k_transfer_net[:,level], 0.0)

            k_transfer_cond_smooth[:,level] = nan_to_value(k_transfer_cond_smooth[:,level], 0.0)
            k_transfer_cond_net_smooth[:,level] = nan_to_value(k_transfer_cond_net_smooth[:,level], 0.0)
            k_transfer_smooth[:,level] = nan_to_value(k_transfer_smooth[:,level], 0.0)
            k_transfer_net_smooth[:,level] = nan_to_value(k_transfer_net_smooth[:,level], 0.0)

            tau_transfer_cond[:,level] = 1 / k_transfer_cond[:,level]
            tau_transfer_cond_net[:,level] = 1 / k_transfer_cond_net[:,level]
            tau_transfer[:,level] = 1 / k_transfer[:,level]
            tau_transfer_net[:,level] = 1 / k_transfer_net[:,level]

            tau_transfer_cond_smooth[:,level] = 1 / k_transfer_cond_smooth[:,level]
            tau_transfer_cond_net_smooth[:,level] = 1 / k_transfer_cond_net_smooth[:,level]
            tau_transfer_smooth[:,level] = 1 / k_transfer_smooth[:,level]
            tau_transfer_net_smooth[:,level] = 1 / k_transfer_net_smooth[:,level]

            big_value = 1e50

            tau_transfer_cond[:,level] = inf_to_value(tau_transfer_cond[:,level], big_value)
            tau_transfer_cond_net[:,level] = inf_to_value(tau_transfer_cond_net[:,level], big_value)
            tau_transfer[:,level] = inf_to_value(tau_transfer[:,level], big_value)
            tau_transfer_net[:,level] = inf_to_value(tau_transfer_net[:,level], big_value)

            tau_transfer_cond_smooth[:,level] = inf_to_value(tau_transfer_cond_smooth[:,level], big_value)
            tau_transfer_cond_net_smooth[:,level] = inf_to_value(tau_transfer_cond_net_smooth[:,level], big_value)
            tau_transfer_smooth[:,level] = inf_to_value(tau_transfer_smooth[:,level], big_value)
            tau_transfer_net_smooth[:,level] = inf_to_value(tau_transfer_net_smooth[:,level], big_value)

            tau_day[level] = mean_day(time[1:], tau_transfer[:,level])
            tau_day_cond[level] = mean_day(time[1:], tau_transfer_cond[:,level])
            tau_night[level] = mean_night(time[1:], tau_transfer[:,level])
            tau_night_cond[level] = mean_night(time[1:], tau_transfer_cond[:,level])

        print "Writing data %s %s..." % (coag_suffix, type_suffix)
        filename = "%s/aging_%s_%s_%%s.txt" % (data_prefix, coag_suffix, type_suffix)

        savetxt(filename % "a_smooth", data_a_smooth, fmt = "%.20e")
        savetxt(filename % "f_smooth", data_f_smooth, fmt = "%.20e")
        savetxt(filename % "emit_a_smooth", data_emit_a_smooth, fmt = "%.20e")
        savetxt(filename % "emit_f_smooth", data_emit_f_smooth, fmt = "%.20e")
        savetxt(filename % "dilution_a_smooth", data_dilution_a_smooth, fmt = "%.20e")
        savetxt(filename % "dilution_f_smooth", data_dilution_f_smooth, fmt = "%.20e")
        savetxt(filename % "halving_a_smooth", data_halving_a_smooth, fmt = "%.20e")
        savetxt(filename % "halving_f_smooth", data_halving_f_smooth, fmt = "%.20e")
        savetxt(filename % "cond_a_a_smooth", data_cond_a_a_smooth, fmt = "%.20e")
        savetxt(filename % "cond_a_f_smooth", data_cond_a_f_smooth, fmt = "%.20e")
        savetxt(filename % "cond_f_a_smooth", data_cond_f_a_smooth, fmt = "%.20e")
        savetxt(filename % "cond_f_f_smooth", data_cond_f_f_smooth, fmt = "%.20e")
        savetxt(filename % "coag_gain_a_smooth", data_coag_gain_a_smooth, fmt = "%.20e")
        savetxt(filename % "coag_gain_f_smooth", data_coag_gain_f_smooth, fmt = "%.20e")
        savetxt(filename % "coag_loss_a_a_smooth", data_coag_loss_a_a_smooth, fmt = "%.20e")
        savetxt(filename % "coag_loss_a_f_smooth", data_coag_loss_a_f_smooth, fmt = "%.20e")
        savetxt(filename % "coag_loss_f_a_smooth", data_coag_loss_f_a_smooth, fmt = "%.20e")
        savetxt(filename % "coag_loss_f_f_smooth", data_coag_loss_f_f_smooth, fmt = "%.20e")
        savetxt(filename % "k_transfer_cond", k_transfer_cond, fmt = "%.20e")
        savetxt(filename % "k_transfer_cond_net", k_transfer_cond_net, fmt = "%.20e")
        savetxt(filename % "k_transfer", k_transfer, fmt = "%.20e")
        savetxt(filename % "k_transfer_net", k_transfer_net, fmt = "%.20e")
        savetxt(filename % "k_transfer_cond_smooth", k_transfer_cond_smooth, fmt = "%.20e")
        savetxt(filename % "k_transfer_cond_net_smooth", k_transfer_cond_net_smooth, fmt = "%.20e")
        savetxt(filename % "k_transfer_smooth", k_transfer_smooth, fmt = "%.20e")
        savetxt(filename % "k_transfer_net_smooth", k_transfer_net_smooth, fmt = "%.20e")
        savetxt(filename % "tau_transfer_cond", tau_transfer_cond, fmt = "%.20e")
        savetxt(filename % "tau_transfer_cond_net", tau_transfer_cond_net, fmt = "%.20e")
        savetxt(filename % "tau_transfer", tau_transfer, fmt = "%.20e")
        savetxt(filename % "tau_transfer_net", tau_transfer_net, fmt = "%.20e")
        savetxt(filename % "tau_transfer_cond_smooth", tau_transfer_cond_smooth, fmt = "%.20e")
        savetxt(filename % "tau_transfer_cond_net_smooth", tau_transfer_cond_net_smooth, fmt = "%.20e")
        savetxt(filename % "tau_transfer_smooth", tau_transfer_smooth, fmt = "%.20e")
        savetxt(filename % "tau_transfer_net_smooth", tau_transfer_net_smooth, fmt = "%.20e")
        savetxt(filename % "tau_day", tau_day, fmt = "%.20e")
        savetxt(filename % "tau_day_cond", tau_day_cond, fmt = "%.20e")
        savetxt(filename % "tau_night", tau_night, fmt = "%.20e")
        savetxt(filename % "tau_night_cond", tau_night_cond, fmt = "%.20e")

print "max_error_num: ", max_error_num
print "max_error_mass: ", max_error_mass
