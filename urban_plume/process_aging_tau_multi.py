#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, pyx
sys.path.append("../tool")
from pmc_pyx import *
sys.path.append(".")
from process_aging_config import *
from numpy import *

def delta(arr):
    return (arr[1:] - arr[:-1])

def mid(arr):
    return (arr[1:] + arr[:-1]) / 2.0

def filter_inf(plot_data):
    return [[x, y] for [x, y] in plot_data if isfinite(y)]

def sign(x):
    if x > 0:
        return 1
    elif x < 0:
        return -1
    return 0

def mean(x):
    return float(sum(x)) / len(x)

def mean_day(time, data):
    return mean([x for [t,x] in filter_inf(zip(time, data))
                 if t >= 6 * 3600 - 0.1 and t <= 9 * 3600 + 0.1])

def mean_night(time, data):
    return mean([x for [t,x] in filter_inf(zip(time, data))
                 if t >= 12 * 3600 - 0.1 and t <= 22 * 3600 + 0.1])

def chop_sign_data_helper(plot_data, chopped_data):
    if len(plot_data) == 0:
        return chopped_data
    if sign(plot_data[0][1]) == 0:
        return chop_sign_data_helper(plot_data[1:], chopped_data)
    i = 0
    while (i < len(plot_data)) \
              and (sign(plot_data[i][1]) == sign(plot_data[0][1])):
        i += 1
    chopped_data.append(plot_data[0:i])
    return chop_sign_data_helper(plot_data[i:], chopped_data)

def chop_sign_data(plot_data):
    return chop_sign_data_helper(plot_data, [])

def smooth(x,window_len=10,window='hanning'):
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len<3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is not one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    if window == 'flat':
        w=ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

dilution_rate = 1.5e-5 # s^{-1}

data_prefix = "out"

smooth_window_len = 60

grey_level = 0.2

max_error_num = 0
max_error_mass = 0.0

def process_data(level, output_means = False, output_plots = False, outf = None):
    global max_error_num, max_error_mass
    for coag in [True, False]:
        for num in [True, False]:
            if coag:
                coag_suffix = "wc"
            else:
                coag_suffix = "nc"

            if num:
                type_suffix = "num"
                type = int
            else:
                type_suffix = "mass"
                type = float

            filename_pattern = "aging_%s_([0-9]{8})_time.txt" % coag_suffix
            filename_re = re.compile(filename_pattern)
            filename_list = get_filename_list(data_prefix, filename_pattern)
            n = len(filename_list)

            time = zeros(n, float)
            height = zeros(n, float)
            comp_vol = zeros(n, float)

            data_a = zeros(n, type)
            data_f = zeros(n, type)
            data_emit_a = zeros(n - 1, type)
            data_emit_f = zeros(n - 1, type)
            data_dilution_a = zeros(n - 1, type)
            data_dilution_f = zeros(n - 1, type)
            data_halving_a = zeros(n - 1, type)
            data_halving_f = zeros(n - 1, type)
            data_cond_a_a = zeros(n - 1, type)
            data_cond_a_f = zeros(n - 1, type)
            data_cond_f_a = zeros(n - 1, type)
            data_cond_f_f = zeros(n - 1, type)
            data_coag_gain_a = zeros(n - 1, type)
            data_coag_gain_f = zeros(n - 1, type)
            data_coag_loss_a_a = zeros(n - 1, type)
            data_coag_loss_a_f = zeros(n - 1, type)
            data_coag_loss_f_a = zeros(n - 1, type)
            data_coag_loss_f_f = zeros(n - 1, type)

            first_time = True
            for (i, filename) in enumerate(filename_list):
                match = filename_re.search(filename)
                if not match:
                    raise Exception()
                key = match.group(1)

                filename_base = "%s/aging_%s_%s" % (data_prefix, coag_suffix, key)
                filename_base_type = "%s/aging_%s_%s_%s" % (data_prefix, coag_suffix, key, type_suffix)

                time_array = loadtxt("%s_time.txt" % filename_base)
                height_array = loadtxt("%s_height.txt" % filename_base)
                comp_vol_array = loadtxt("%s_comp_vol.txt" % filename_base)

                data_array = loadtxt("%s.txt" % filename_base_type, type)
                data_emit_array = loadtxt("%s_emit.txt" % filename_base_type, type)
                data_dilution_array = loadtxt("%s_dilution.txt" % filename_base_type, type)
                data_halving_array = loadtxt("%s_halving.txt" % filename_base_type, type)
                data_cond_array = loadtxt("%s_cond.txt" % filename_base_type, type)
                data_coag_gain_array = loadtxt("%s_coag_gain.txt" % filename_base_type, type)
                data_coag_loss_array = loadtxt("%s_coag_loss.txt" % filename_base_type, type)

                time[i] = float(time_array)
                height[i] = float(height_array)
                comp_vol[i] = float(comp_vol_array)

                k = level + 1
                data_a[i] = data_array[:k].sum()
                data_f[i] = data_array[k:].sum()
                if not first_time:
                    data_emit_a[i-1] = data_emit_array[:k].sum()
                    data_emit_f[i-1] = data_emit_array[k:].sum()
                    data_dilution_a[i-1] = data_dilution_array[:k].sum()
                    data_dilution_f[i-1] = data_dilution_array[k:].sum()
                    data_halving_a[i-1] = data_halving_array[:k].sum()
                    data_halving_f[i-1] = data_halving_array[k:].sum()
                    data_cond_a_a[i-1] = data_cond_array[:k,:k].sum().sum()
                    data_cond_a_f[i-1] = data_cond_array[:k,k:].sum().sum()
                    data_cond_f_a[i-1] = data_cond_array[k:,:k].sum().sum()
                    data_cond_f_f[i-1] = data_cond_array[k:,k:].sum().sum()
                    data_coag_gain_a[i-1] = data_coag_gain_array[:k].sum()
                    data_coag_gain_f[i-1] = data_coag_gain_array[k:].sum()
                    data_coag_loss_a_a[i-1] = data_coag_loss_array[:k,:k].sum().sum()
                    data_coag_loss_a_f[i-1] = data_coag_loss_array[:k,k:].sum().sum()
                    data_coag_loss_f_a[i-1] = data_coag_loss_array[k:,:k].sum().sum()
                    data_coag_loss_f_f[i-1] = data_coag_loss_array[k:,k:].sum().sum()

                first_time = False

            #print "data_a: ", data_a
            #print "data_f: ", data_f
            #print "data_emit_a: ", data_emit_a
            #print "data_emit_f: ", data_emit_f
            #print "data_dilution_a: ", data_dilution_a
            #print "data_dilution_f: ", data_dilution_f
            #print "data_halving_a: ", data_halving_a
            #print "data_halving_f: ", data_halving_f
            #print "data_cond_a_a: ", data_cond_a_a
            #print "data_cond_a_f: ", data_cond_a_f
            #print "data_cond_f_a: ", data_cond_f_a
            #print "data_cond_f_f: ", data_cond_f_f
            #print "data_coag_gain_a: ", data_coag_gain_a
            #print "data_coag_gain_f: ", data_coag_gain_f
            #print "data_coag_loss_a_a: ", data_coag_loss_a_a
            #print "data_coag_loss_a_f: ", data_coag_loss_a_f
            #print "data_coag_loss_f_a: ", data_coag_loss_f_a
            #print "data_coag_loss_f_f: ", data_coag_loss_f_f
            
            max_time_min = max(time) / 60
            start_time_of_day_min = 6 * 60

            data_a_smooth = smooth(data_a, window_len = smooth_window_len)
            data_f_smooth = smooth(data_f, window_len = smooth_window_len)
            data_emit_a_smooth = smooth(data_emit_a, window_len = smooth_window_len)
            data_emit_f_smooth = smooth(data_emit_f, window_len = smooth_window_len)
            data_dilution_a_smooth = smooth(data_dilution_a, window_len = smooth_window_len)
            data_dilution_f_smooth = smooth(data_dilution_f, window_len = smooth_window_len)
            data_halving_a_smooth = smooth(data_halving_a, window_len = smooth_window_len)
            data_halving_f_smooth = smooth(data_halving_f, window_len = smooth_window_len)
            data_cond_a_a_smooth = smooth(data_cond_a_a, window_len = smooth_window_len)
            data_cond_a_f_smooth = smooth(data_cond_a_f, window_len = smooth_window_len)
            data_cond_f_a_smooth = smooth(data_cond_f_a, window_len = smooth_window_len)
            data_cond_f_f_smooth = smooth(data_cond_f_f, window_len = smooth_window_len)
            data_coag_gain_a_smooth = smooth(data_coag_gain_a, window_len = smooth_window_len)
            data_coag_gain_f_smooth = smooth(data_coag_gain_f, window_len = smooth_window_len)
            data_coag_loss_a_a_smooth = smooth(data_coag_loss_a_a, window_len = smooth_window_len)
            data_coag_loss_a_f_smooth = smooth(data_coag_loss_a_f, window_len = smooth_window_len)
            data_coag_loss_f_a_smooth = smooth(data_coag_loss_f_a, window_len = smooth_window_len)
            data_coag_loss_f_f_smooth = smooth(data_coag_loss_f_f, window_len = smooth_window_len)

            delta_data_a = delta(data_a)
            delta_data_f = delta(data_f)

            aged_kp1 = data_a[1:] - (data_cond_a_a + data_cond_f_a +
                                    data_emit_a + data_coag_gain_a)
            aged_k = data_a[:-1] - (data_cond_a_a + data_cond_a_f +
                                   data_dilution_a + data_halving_a +
                                   data_coag_loss_a_a + data_coag_loss_a_f)
            aged_delta = delta_data_a - (data_emit_a - data_dilution_a -
                                        data_halving_a + data_cond_f_a -
                                        data_cond_a_f + data_coag_gain_a -
                                        data_coag_loss_a_a - data_coag_loss_a_f)
            fresh_kp1 = data_f[1:] - (data_cond_a_f + data_cond_f_f +
                                     data_emit_f + data_coag_gain_f)
            fresh_k = data_f[:-1] - (data_cond_f_a + data_cond_f_f +
                                    data_dilution_f + data_halving_f +
                                    data_coag_loss_f_a + data_coag_loss_f_f)
            fresh_delta = delta_data_f - (data_emit_f - data_dilution_f -
                                         data_halving_f + data_cond_a_f -
                                         data_cond_f_a + data_coag_gain_f -
                                         data_coag_loss_f_a - data_coag_loss_f_f)

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
            if num:
                max_error_num = max(max_error_num, max_error)
            else:
                max_error_mass = max(max_error_mass, max_error)

            k_transfer_cond = data_cond_f_a / delta(time) / data_f[1:]
            k_transfer_cond_net = (data_cond_f_a - data_cond_a_f) / delta(time) / data_f[1:]
            k_transfer = (data_cond_f_a + data_coag_loss_f_a) / delta(time) / data_f[1:]
            k_transfer_net = (data_cond_f_a - data_cond_a_f + data_coag_loss_f_a - data_coag_loss_a_f) / delta(time) / data_f[1:]

            k_transfer_cond_smooth = data_cond_f_a_smooth / delta(time) / data_f_smooth[1:]
            k_transfer_cond_net_smooth = (data_cond_f_a_smooth - data_cond_a_f_smooth) / delta(time) / data_f_smooth[1:]
            k_transfer_smooth = (data_cond_f_a_smooth + data_coag_loss_f_a_smooth) / delta(time) / data_f_smooth[1:]
            k_transfer_net_smooth = (data_cond_f_a_smooth - data_cond_a_f_smooth + data_coag_loss_f_a_smooth - data_coag_loss_a_f_smooth) / delta(time) / data_f_smooth[1:]

            tau_transfer_cond = 1 / k_transfer_cond
            tau_transfer_cond_net = 1 / k_transfer_cond_net
            tau_transfer = 1 / k_transfer
            tau_transfer_net = 1 / k_transfer_net

            tau_transfer_cond_smooth = 1 / k_transfer_cond_smooth
            tau_transfer_cond_net_smooth = 1 / k_transfer_cond_net_smooth
            tau_transfer_smooth = 1 / k_transfer_smooth
            tau_transfer_net_smooth = 1 / k_transfer_net_smooth

            if output_means:
                tau_day = mean_day(time[1:], tau_transfer)
                tau_day_cond = mean_day(time[1:], tau_transfer_cond)
                tau_night = mean_night(time[1:], tau_transfer)
                tau_night_cond = mean_night(time[1:], tau_transfer_cond)
                outf.write("%.20e %.20e %.20e %.20e " % (tau_day, tau_day_cond, tau_night, tau_night_cond))

            if not output_plots:
                continue

            ######################################################################

            g = pyx.graph.graphxy(
                width = 10,
                x = pyx.graph.axis.linear(title = r"time (s)",
                                      painter = grid_painter),
                y = pyx.graph.axis.linear(painter = grid_painter),
                key = pyx.graph.key.key(pos = "tr"))

            g.plot(
                pyx.graph.data.points(zip(time, data_a), x = 1, y = 2,
                                      title = "%s a" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

            g.plot(
                pyx.graph.data.points(zip(time, data_f), x = 1, y = 2,
                                      title = "%s f" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

            g.writePDFfile("%s/aging_%s_%s_%d_total.pdf" % (data_prefix, coag_suffix, type_suffix, level))

            ######################################################################

            g = pyx.graph.graphxy(
                width = 10,
                x = pyx.graph.axis.linear(title = r"time (s)",
                                      painter = grid_painter),
                y = pyx.graph.axis.linear(painter = grid_painter),
                key = pyx.graph.key.key(pos = "tr"))

            g.plot(
                pyx.graph.data.points(zip(time[1:], data_emit_a), x = 1, y = 2,
                                      title = "%s emit a" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

            g.plot(
                pyx.graph.data.points(zip(time[1:], data_emit_f), x = 1, y = 2,
                                      title = "%s emit f" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

            g.writePDFfile("%s/aging_%s_%s_%d_emit.pdf" % (data_prefix, coag_suffix, type_suffix, level))

            ######################################################################

            g = pyx.graph.graphxy(
                width = 10,
                x = pyx.graph.axis.linear(title = r"time (s)",
                                      painter = grid_painter),
                y = pyx.graph.axis.linear(painter = grid_painter),
                key = pyx.graph.key.key(pos = "tr"))

            g.plot(
                pyx.graph.data.points(zip(time[1:], data_dilution_a), x = 1, y = 2,
                                      title = "%s dilution a" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

            g.plot(
                pyx.graph.data.points(zip(time[1:], data_dilution_f), x = 1, y = 2,
                                      title = "%s dilution f" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

            g.writePDFfile("%s/aging_%s_%s_%d_dilution.pdf" % (data_prefix, coag_suffix, type_suffix, level))

            ######################################################################

            g = pyx.graph.graphxy(
                width = 10,
                x = pyx.graph.axis.linear(title = r"time (s)",
                                      painter = grid_painter),
                y = pyx.graph.axis.linear(painter = grid_painter),
                key = pyx.graph.key.key(pos = "tr"))

            g.plot(
                pyx.graph.data.points(zip(time[1:], data_halving_a), x = 1, y = 2,
                                      title = "%s halving a" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

            g.plot(
                pyx.graph.data.points(zip(time[1:], data_halving_f), x = 1, y = 2,
                                      title = "%s halving f" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

            g.writePDFfile("%s/aging_%s_%s_%d_halving.pdf" % (data_prefix, coag_suffix, type_suffix, level))

            ######################################################################

            g = pyx.graph.graphxy(
                width = 10,
                x = pyx.graph.axis.linear(title = r"time (s)",
                                      painter = grid_painter),
                y = pyx.graph.axis.linear(painter = grid_painter),
                key = pyx.graph.key.key(pos = "tr"))

            g.plot(
                pyx.graph.data.points(zip(time[1:], data_cond_a_a), x = 1, y = 2,
                                      title = "%s cond a a" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

            g.plot(
                pyx.graph.data.points(zip(time[1:], data_cond_f_f), x = 1, y = 2,
                                      title = "%s cond f f" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

            g.writePDFfile("%s/aging_%s_%s_%d_cond_remain.pdf" % (data_prefix, coag_suffix, type_suffix, level))

            ######################################################################

            g = pyx.graph.graphxy(
                width = 10,
                x = pyx.graph.axis.linear(title = r"time (s)",
                                      painter = grid_painter),
                y = pyx.graph.axis.linear(painter = grid_painter),
                key = pyx.graph.key.key(pos = "tr"))

            g.plot(
                pyx.graph.data.points(zip(time[1:], data_cond_a_f), x = 1, y = 2,
                                      title = "%s cond a f" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

            g.plot(
                pyx.graph.data.points(zip(time[1:], data_cond_f_a), x = 1, y = 2,
                                      title = "%s cond f a" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

            g.writePDFfile("%s/aging_%s_%s_%d_cond_transfer.pdf" % (data_prefix, coag_suffix, type_suffix, level))

            ######################################################################

            if max(data_coag_gain_a) > 0 or max(data_coag_gain_f) > 0:
                g = pyx.graph.graphxy(
                    width = 10,
                    x = pyx.graph.axis.linear(title = r"time (s)",
                                          painter = grid_painter),
                    y = pyx.graph.axis.linear(painter = grid_painter),
                    key = pyx.graph.key.key(pos = "tr"))

                g.plot(
                    pyx.graph.data.points(zip(time[1:], data_coag_gain_a), x = 1, y = 2,
                                          title = "%s coag gain a" % type_suffix),
                    styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

                g.plot(
                    pyx.graph.data.points(zip(time[1:], data_coag_gain_f), x = 1, y = 2,
                                          title = "%s coag gain f" % type_suffix),
                    styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

                g.writePDFfile("%s/aging_%s_%s_%d_coag_gain.pdf" % (data_prefix, coag_suffix, type_suffix, level))

            ######################################################################

            if max(data_coag_loss_a_a) > 0 \
                   or max(data_coag_loss_f_f) > 0:
                g = pyx.graph.graphxy(
                    width = 10,
                    x = pyx.graph.axis.linear(title = r"time (s)",
                                          painter = grid_painter),
                    y = pyx.graph.axis.linear(painter = grid_painter),
                    key = pyx.graph.key.key(pos = "tr"))

                g.plot(
                    pyx.graph.data.points(zip(time[1:], data_coag_loss_a_a), x = 1, y = 2,
                                          title = "%s coag loss a a" % type_suffix),
                    styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

                g.plot(
                    pyx.graph.data.points(zip(time[1:], data_coag_loss_f_f), x = 1, y = 2,
                                          title = "%s coag loss f f" % type_suffix),
                    styles = [pyx.graph.style.line(lineattrs = [color_list[3]])])

                g.writePDFfile("%s/aging_%s_%s_%d_coag_loss_remain.pdf" % (data_prefix, coag_suffix, type_suffix, level))

            ######################################################################

            if max(data_coag_loss_a_f) > 0 \
                   or max(data_coag_loss_f_a) > 0:
                g = pyx.graph.graphxy(
                    width = 10,
                    x = pyx.graph.axis.linear(title = r"time (s)",
                                          painter = grid_painter),
                    y = pyx.graph.axis.linear(painter = grid_painter),
                    key = pyx.graph.key.key(pos = "tr"))

                g.plot(
                    pyx.graph.data.points(zip(time[1:], data_coag_loss_a_f), x = 1, y = 2,
                                          title = "%s coag loss a f" % type_suffix),
                    styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

                g.plot(
                    pyx.graph.data.points(zip(time[1:], data_coag_loss_f_a), x = 1, y = 2,
                                          title = "%s coag loss f a" % type_suffix),
                    styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

                g.writePDFfile("%s/aging_%s_%s_%d_coag_loss_transfer.pdf" % (data_prefix, coag_suffix, type_suffix, level))

            ######################################################################

            g = pyx.graph.graphxy(
                width = 10,
                x = pyx.graph.axis.linear(title = r"time (s)",
                                      painter = grid_painter),
                y = pyx.graph.axis.linear(painter = grid_painter),
                key = pyx.graph.key.key(pos = "tr"))

            g.plot(
                pyx.graph.data.points(zip(time, height), x = 1, y = 2,
                                      title = "height"),
                styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

            g.writePDFfile("%s/aging_%s_%s_%d_height.pdf" % (data_prefix, coag_suffix, type_suffix, level))

            ######################################################################

            g = pyx.graph.graphxy(
                width = 10,
                x = pyx.graph.axis.linear(title = r"time (s)",
                                      painter = grid_painter),
                y = pyx.graph.axis.linear(painter = grid_painter),
                key = pyx.graph.key.key(pos = "tr"))

            g.plot(
                pyx.graph.data.points(zip(time, comp_vol), x = 1, y = 2,
                                      title = "comp vol"),
                styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

            g.writePDFfile("%s/aging_%s_%s_%d_comp_vol.pdf" % (data_prefix, coag_suffix, type_suffix, level))

            ######################################################################

            g = pyx.graph.graphxy(
                width = 10,
                x = pyx.graph.axis.linear(title = r"time (s)",
                                      painter = grid_painter),
                y = pyx.graph.axis.linear(painter = grid_painter),
                key = pyx.graph.key.key(pos = "tr"))

            g.plot(
                pyx.graph.data.points(zip(time[1:], k_transfer), x = 1, y = 2,
                                      title = "%s k transfer" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

            g.plot(
                pyx.graph.data.points(zip(time[1:], k_transfer_net), x = 1, y = 2,
                                      title = "%s k transfer net" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

            g.plot(
                pyx.graph.data.points(zip(time[1:], k_transfer_cond), x = 1, y = 2,
                                      title = "%s k transfer cond" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[2]])])

            g.plot(
                pyx.graph.data.points(zip(time[1:], k_transfer_cond_net), x = 1, y = 2,
                                      title = "%s k transfer cond net" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[3]])])

            g.writePDFfile("%s/aging_%s_%s_%d_k.pdf" % (data_prefix, coag_suffix, type_suffix, level))

            ######################################################################

            g = pyx.graph.graphxy(
                width = 10,
                x = pyx.graph.axis.linear(title = r"time (s)",
                                      painter = grid_painter),
                y = pyx.graph.axis.linear(painter = grid_painter),
                key = pyx.graph.key.key(pos = "tr"))

            g.plot(
                pyx.graph.data.points(filter_inf(zip(time[1:], tau_transfer)), x = 1, y = 2,
                                      title = "%s tau transfer" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

            g.plot(
                pyx.graph.data.points(filter_inf(zip(time[1:], tau_transfer_net)), x = 1, y = 2,
                                      title = "%s tau transfer net" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

            g.plot(
                pyx.graph.data.points(filter_inf(zip(time[1:], tau_transfer_cond)), x = 1, y = 2,
                                      title = "%s tau transfer cond" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[2]])])

            g.plot(
                pyx.graph.data.points(filter_inf(zip(time[1:], tau_transfer_cond_net)), x = 1, y = 2,
                                      title = "%s tau transfer cond net" % type_suffix),
                styles = [pyx.graph.style.line(lineattrs = [color_list[3]])])

            g.writePDFfile("%s/aging_%s_%s_%d_tau.pdf" % (data_prefix, coag_suffix, type_suffix, level))

            ######################################################################

            g = pyx.graph.graphxy(
                width = 10,
                x = pyx.graph.axis.linear(min = 0.,
                                          max = max_time_min,
                                          parter = graph.axis.parter.linear(tickdists
                                                                            = [6 * 60, 3 * 60]),
                                          texter = time_of_day(base_time
                                                               = start_time_of_day_min),
                                          title = "local standard time (LST) (hours:minutes)",
                                          painter = grid_painter),
                y = pyx.graph.axis.linear(min = -2e-4,
                                          max = 4e-3,
                                          title = r"aging rate $k$ ($\rm s^{-1}$)",
                                          painter = grid_painter),
                key = pyx.graph.key.key(pos = "tr"))

            for (i, k_data) in enumerate([k_transfer, k_transfer_net,
                                          k_transfer_cond, k_transfer_cond_net]):
                plot_data = filter_inf(zip(time[1:] / 60, k_data))
                grey_color = pyx.color.hsb(color_list[i].hsb().color["h"], grey_level, 1)
                g.plot(
                    pyx.graph.data.points(plot_data, x = 1, y = 2, title = None),
                    styles = [pyx.graph.style.line(lineattrs = [grey_color])])

            for (i, (k_data, title)) in enumerate(zip([k_transfer_smooth, k_transfer_net_smooth,
                                                       k_transfer_cond_smooth, k_transfer_cond_net_smooth],
                                                      ["%s k transfer smooth" % type_suffix, "%s k transfer net smooth" % type_suffix,
                                                       "%s k transfer cond smooth" % type_suffix, "%s k transfer cond net smooth" % type_suffix])):
                plot_data = filter_inf(zip(time[1:] / 60, k_data))
                g.plot(
                    pyx.graph.data.points(plot_data, x = 1, y = 2, title = title),
                    styles = [pyx.graph.style.line(lineattrs = [color_list[i]])])

            g.writePDFfile("%s/aging_%s_%s_%d_k.pdf" % (data_prefix, coag_suffix, type_suffix, level))

            ######################################################################

            g = pyx.graph.graphxy(
                width = 10,
                x = pyx.graph.axis.linear(min = 0.,
                                          max = max_time_min,
                                          parter = graph.axis.parter.linear(tickdists
                                                                            = [6 * 60, 3 * 60]),
                                          texter = time_of_day(base_time
                                                               = start_time_of_day_min),
                                          title = "local standard time (LST) (hours:minutes)",
                                          painter = grid_painter),
                y = pyx.graph.axis.log(min = 1e-2,
                                       max = 1e3,
                                       title = r"aging timescale $\tau$ (hours)",
                                       painter = grid_painter))

            for (i, tau_data) in enumerate([tau_transfer, tau_transfer_net,
                                            tau_transfer_cond, tau_transfer_cond_net]):
                plot_data = filter_inf(zip(time[1:] / 60, tau_data / 3600))
                chopped_data = chop_sign_data(plot_data)
                grey_color = pyx.color.hsb(color_list[i].hsb().color["h"], grey_level, 1)
                for signed_data in chopped_data:
                    if signed_data[0][1] > 0:
                        g.plot(
                            pyx.graph.data.points(signed_data, x = 1, y = 2),
                            styles = [pyx.graph.style.line(lineattrs = [grey_color])])
                    else:
                        signed_data = [[t,-d] for [t,d] in signed_data]
                        g.plot(
                            pyx.graph.data.points(signed_data, x = 1, y = 2),
                            styles = [pyx.graph.style.line(lineattrs = [grey_color,
                                                                        pyx.style.linestyle.dashed])])

            for (i, tau_data) in enumerate([tau_transfer_smooth, tau_transfer_net_smooth,
                                            tau_transfer_cond_smooth, tau_transfer_cond_net_smooth]):
                plot_data = filter_inf(zip(time[1:] / 60, tau_data / 3600))
                chopped_data = chop_sign_data(plot_data)
                for signed_data in chopped_data:
                    if signed_data[0][1] > 0:
                        g.plot(
                            pyx.graph.data.points(signed_data, x = 1, y = 2),
                            styles = [pyx.graph.style.line(lineattrs = [color_list[i]])])
                    else:
                        signed_data = [[t,-d] for [t,d] in signed_data]
                        g.plot(
                            pyx.graph.data.points(signed_data, x = 1, y = 2),
                            styles = [pyx.graph.style.line(lineattrs = [color_list[i],
                                                                        pyx.style.linestyle.dashed])])

            g.writePDFfile("%s/aging_%s_%s_%d_tau_log.pdf" % (data_prefix, coag_suffix, type_suffix, level))

for plot_ss in [0.001, 0.003, 0.006, 0.01]:
    level = int(ss_active_axis.closest_edge(array(plot_ss)))
    print "plot_ss = ", plot_ss
    print "level = ", level
    process_data(level, output_means = False, output_plots = True)

outf_means = open("%s/aging_tau_means.txt" % data_prefix, "w")
for level in range(ss_active_axis.n_bin + 1):
    print "plot_ss = ", ss_active_axis.edge(level)
    print "level = ", level
    outf_means.write("%d %.20e " % (level, ss_active_axis.edge(level)))
    process_data(level, output_means = True, output_plots = False, outf = outf_means)
    outf_means.write("\n")
outf_means.close()

print "max_error_num: ", max_error_num
print "max_error_mass: ", max_error_mass
