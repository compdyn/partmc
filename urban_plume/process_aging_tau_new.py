#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, pyx
from numpy import *
sys.path.append("../tool")
from pmc_pyx import *

def delta(arr):
    return (arr[1:] - arr[:-1])

def mid(arr):
    return (arr[1:] + arr[:-1]) / 2.0

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
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

dilution_rate = 1.5e-5 # s^{-1}

data_prefix = "aging_data/1"

smooth_window_len = 60

level = 4

grey_level = 0.8

for coag in [True, False]:
    if coag:
        coag_suffix = "wc"
    else:
        coag_suffix = "nc"

    height_array = loadtxt("%s/aging_%s_height.txt" % (data_prefix, coag_suffix), unpack = True)
    comp_vol_array = loadtxt("%s/aging_%s_comp_vol.txt" % (data_prefix, coag_suffix), unpack = True)

    num_a_array = loadtxt("%s/aging_%s_num_a.txt" % (data_prefix, coag_suffix), unpack = True)
    num_f_array = loadtxt("%s/aging_%s_num_f.txt" % (data_prefix, coag_suffix), unpack = True)
    num_emit_a_array = loadtxt("%s/aging_%s_num_emit_a.txt" % (data_prefix, coag_suffix), unpack = True)
    num_emit_f_array = loadtxt("%s/aging_%s_num_emit_f.txt" % (data_prefix, coag_suffix), unpack = True)
    num_dilution_a_array = loadtxt("%s/aging_%s_num_dilution_a.txt" % (data_prefix, coag_suffix), unpack = True)
    num_dilution_f_array = loadtxt("%s/aging_%s_num_dilution_f.txt" % (data_prefix, coag_suffix), unpack = True)
    num_halving_a_array = loadtxt("%s/aging_%s_num_halving_a.txt" % (data_prefix, coag_suffix), unpack = True)
    num_halving_f_array = loadtxt("%s/aging_%s_num_halving_f.txt" % (data_prefix, coag_suffix), unpack = True)
    num_cond_a_a_array = loadtxt("%s/aging_%s_num_cond_a_a.txt" % (data_prefix, coag_suffix), unpack = True)
    num_cond_a_f_array = loadtxt("%s/aging_%s_num_cond_a_f.txt" % (data_prefix, coag_suffix), unpack = True)
    num_cond_f_a_array = loadtxt("%s/aging_%s_num_cond_f_a.txt" % (data_prefix, coag_suffix), unpack = True)
    num_cond_f_f_array = loadtxt("%s/aging_%s_num_cond_f_f.txt" % (data_prefix, coag_suffix), unpack = True)
    num_coag_gain_a_array = loadtxt("%s/aging_%s_num_coag_gain_a.txt" % (data_prefix, coag_suffix), unpack = True)
    num_coag_gain_f_array = loadtxt("%s/aging_%s_num_coag_gain_f.txt" % (data_prefix, coag_suffix), unpack = True)
    num_coag_loss_a_a_array = loadtxt("%s/aging_%s_num_coag_loss_a_a.txt" % (data_prefix, coag_suffix), unpack = True)
    num_coag_loss_a_f_array = loadtxt("%s/aging_%s_num_coag_loss_a_f.txt" % (data_prefix, coag_suffix), unpack = True)
    num_coag_loss_f_a_array = loadtxt("%s/aging_%s_num_coag_loss_f_a.txt" % (data_prefix, coag_suffix), unpack = True)
    num_coag_loss_f_f_array = loadtxt("%s/aging_%s_num_coag_loss_f_f.txt" % (data_prefix, coag_suffix), unpack = True)


    time = comp_vol_array[0,:]
    max_time_min = max(time) / 60
    start_time_of_day_min = 6 * 60
    height = height_array[1,:]
    comp_vol = comp_vol_array[1,:]

    num_a = num_a_array[level,:]
    num_f = num_f_array[level,:]
    num_emit_a = num_emit_a_array[level,1:]
    num_emit_f = num_emit_f_array[level,1:]
    num_dilution_a = num_dilution_a_array[level,1:]
    num_dilution_f = num_dilution_f_array[level,1:]
    num_halving_a = num_halving_a_array[level,1:]
    num_halving_f = num_halving_f_array[level,1:]
    num_cond_a_a = num_cond_a_a_array[level,1:]
    num_cond_a_f = num_cond_a_f_array[level,1:]
    num_cond_f_a = num_cond_f_a_array[level,1:]
    num_cond_f_f = num_cond_f_f_array[level,1:]
    num_coag_gain_a = num_coag_gain_a_array[level,1:]
    num_coag_gain_f = num_coag_gain_f_array[level,1:]
    num_coag_loss_a_a = num_coag_loss_a_a_array[level,1:]
    num_coag_loss_a_f = num_coag_loss_a_f_array[level,1:]
    num_coag_loss_f_a = num_coag_loss_f_a_array[level,1:]
    num_coag_loss_f_f = num_coag_loss_f_f_array[level,1:]

    delta_num_a = delta(num_a)
    delta_num_f = delta(num_f)

    aged_kp1 = num_a[1:] - (num_cond_a_a + num_cond_f_a +
                            num_emit_a + num_coag_gain_a)
    aged_k = num_a[:-1] - (num_cond_a_a + num_cond_a_f +
                           num_dilution_a + num_halving_a +
                           num_coag_loss_a_a + num_coag_loss_a_f)
    aged_delta = delta_num_a - (num_emit_a - num_dilution_a -
                                num_halving_a + num_cond_f_a -
                                num_cond_a_f + num_coag_gain_a -
                                num_coag_loss_a_a - num_coag_loss_a_f)
    fresh_kp1 = num_f[1:] - (num_cond_a_f + num_cond_f_f +
                             num_emit_f + num_coag_gain_f)
    fresh_k = num_f[:-1] - (num_cond_f_a + num_cond_f_f +
                            num_dilution_f + num_halving_f +
                            num_coag_loss_f_a + num_coag_loss_f_f)
    fresh_delta = delta_num_f - (num_emit_f - num_dilution_f -
                                 num_halving_f + num_cond_a_f -
                                 num_cond_f_a + num_coag_gain_f -
                                 num_coag_loss_f_a - num_coag_loss_f_f)

    print "aged k+1: ", aged_kp1
    print "aged k: ", aged_k
    print "aged delta: ", aged_delta
    print "fresh k+1: ", fresh_kp1
    print "fresh k: ", fresh_k
    print "fresh delta: ", fresh_delta

    print "max(aged k+1): ", max(abs(aged_kp1))
    print "max(aged k): ", max(abs(aged_k))
    print "max(aged delta): ", max(abs(aged_delta))
    print "max(fresh k+1): ", max(abs(fresh_kp1))
    print "max(fresh k): ", max(abs(fresh_k))
    print "max(fresh delta): ", max(abs(fresh_delta))

    ######################################################################

    g = pyx.graph.graphxy(
        width = 10,
        x = pyx.graph.axis.linear(title = r"time (s)",
                              painter = grid_painter),
        y = pyx.graph.axis.linear(painter = grid_painter),
        key = pyx.graph.key.key(pos = "tr"))

    g.plot(
        pyx.graph.data.points(zip(time, num_a), x = 1, y = 2,
                              title = "num a"),
        styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

    g.plot(
        pyx.graph.data.points(zip(time, num_f), x = 1, y = 2,
                              title = "num f"),
        styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

    g.writePDFfile("%s/aging_num_%d_%s.pdf" % (data_prefix, level, coag_suffix))

    ######################################################################

    g = pyx.graph.graphxy(
        width = 10,
        x = pyx.graph.axis.linear(title = r"time (s)",
                              painter = grid_painter),
        y = pyx.graph.axis.linear(painter = grid_painter),
        key = pyx.graph.key.key(pos = "tr"))

    g.plot(
        pyx.graph.data.points(zip(time[1:], num_emit_a), x = 1, y = 2,
                              title = "num emit a"),
        styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

    g.plot(
        pyx.graph.data.points(zip(time[1:], num_emit_f), x = 1, y = 2,
                              title = "num emit f"),
        styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

    g.writePDFfile("%s/aging_emit_%d_%s.pdf" % (data_prefix, level, coag_suffix))

    ######################################################################

    g = pyx.graph.graphxy(
        width = 10,
        x = pyx.graph.axis.linear(title = r"time (s)",
                              painter = grid_painter),
        y = pyx.graph.axis.linear(painter = grid_painter),
        key = pyx.graph.key.key(pos = "tr"))

    g.plot(
        pyx.graph.data.points(zip(time[1:], num_dilution_a), x = 1, y = 2,
                              title = "num dilution a"),
        styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

    g.plot(
        pyx.graph.data.points(zip(time[1:], num_dilution_f), x = 1, y = 2,
                              title = "num dilution f"),
        styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

    g.writePDFfile("%s/aging_dilution_%d_%s.pdf" % (data_prefix, level, coag_suffix))

    ######################################################################

    g = pyx.graph.graphxy(
        width = 10,
        x = pyx.graph.axis.linear(title = r"time (s)",
                              painter = grid_painter),
        y = pyx.graph.axis.linear(painter = grid_painter),
        key = pyx.graph.key.key(pos = "tr"))

    g.plot(
        pyx.graph.data.points(zip(time[1:], num_halving_a), x = 1, y = 2,
                              title = "num halving a"),
        styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

    g.plot(
        pyx.graph.data.points(zip(time[1:], num_halving_f), x = 1, y = 2,
                              title = "num halving f"),
        styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

    g.writePDFfile("%s/aging_halving_%d_%s.pdf" % (data_prefix, level, coag_suffix))

    ######################################################################

    g = pyx.graph.graphxy(
        width = 10,
        x = pyx.graph.axis.linear(title = r"time (s)",
                              painter = grid_painter),
        y = pyx.graph.axis.linear(painter = grid_painter),
        key = pyx.graph.key.key(pos = "tr"))

    g.plot(
        pyx.graph.data.points(zip(time[1:], num_cond_a_a), x = 1, y = 2,
                              title = "num cond a a"),
        styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

    g.plot(
        pyx.graph.data.points(zip(time[1:], num_cond_f_f), x = 1, y = 2,
                              title = "num cond f f"),
        styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

    g.writePDFfile("%s/aging_cond_remain_%d_%s.pdf" % (data_prefix, level, coag_suffix))

    ######################################################################

    g = pyx.graph.graphxy(
        width = 10,
        x = pyx.graph.axis.linear(title = r"time (s)",
                              painter = grid_painter),
        y = pyx.graph.axis.linear(painter = grid_painter),
        key = pyx.graph.key.key(pos = "tr"))

    g.plot(
        pyx.graph.data.points(zip(time[1:], num_cond_a_f), x = 1, y = 2,
                              title = "num cond a f"),
        styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

    g.plot(
        pyx.graph.data.points(zip(time[1:], num_cond_f_a), x = 1, y = 2,
                              title = "num cond f a"),
        styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

    g.writePDFfile("%s/aging_cond_transfer_%d_%s.pdf" % (data_prefix, level, coag_suffix))

    ######################################################################

    if max(num_coag_gain_a) > 0 or max(num_coag_gain_f) > 0:
        g = pyx.graph.graphxy(
            width = 10,
            x = pyx.graph.axis.linear(title = r"time (s)",
                                  painter = grid_painter),
            y = pyx.graph.axis.linear(painter = grid_painter),
            key = pyx.graph.key.key(pos = "tr"))

        g.plot(
            pyx.graph.data.points(zip(time[1:], num_coag_gain_a), x = 1, y = 2,
                                  title = "num coag gain a"),
            styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

        g.plot(
            pyx.graph.data.points(zip(time[1:], num_coag_gain_f), x = 1, y = 2,
                                  title = "num coag gain f"),
            styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

        g.writePDFfile("%s/aging_coag_gain_%d_%s.pdf" % (data_prefix, level, coag_suffix))

    ######################################################################

    if max(num_coag_loss_a_a) > 0 \
           or max(num_coag_loss_f_f) > 0:
        g = pyx.graph.graphxy(
            width = 10,
            x = pyx.graph.axis.linear(title = r"time (s)",
                                  painter = grid_painter),
            y = pyx.graph.axis.linear(painter = grid_painter),
            key = pyx.graph.key.key(pos = "tr"))

        g.plot(
            pyx.graph.data.points(zip(time[1:], num_coag_loss_a_a), x = 1, y = 2,
                                  title = "num coag loss a a"),
            styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

        g.plot(
            pyx.graph.data.points(zip(time[1:], num_coag_loss_f_f), x = 1, y = 2,
                                  title = "num coag loss f f"),
            styles = [pyx.graph.style.line(lineattrs = [color_list[3]])])

        g.writePDFfile("%s/aging_coag_loss_remain_%d_%s.pdf" % (data_prefix, level, coag_suffix))

    ######################################################################

    if max(num_coag_loss_a_f) > 0 \
           or max(num_coag_loss_f_a) > 0:
        g = pyx.graph.graphxy(
            width = 10,
            x = pyx.graph.axis.linear(title = r"time (s)",
                                  painter = grid_painter),
            y = pyx.graph.axis.linear(painter = grid_painter),
            key = pyx.graph.key.key(pos = "tr"))

        g.plot(
            pyx.graph.data.points(zip(time[1:], num_coag_loss_a_f), x = 1, y = 2,
                                  title = "num coag loss a f"),
            styles = [pyx.graph.style.line(lineattrs = [color_list[0]])])

        g.plot(
            pyx.graph.data.points(zip(time[1:], num_coag_loss_f_a), x = 1, y = 2,
                                  title = "num coag loss f a"),
            styles = [pyx.graph.style.line(lineattrs = [color_list[1]])])

        g.writePDFfile("%s/aging_coag_loss_transfer_%d_%s.pdf" % (data_prefix, level, coag_suffix))

    ######################################################################
