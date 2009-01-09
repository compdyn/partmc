#!/usr/bin/env python
# Copyright (C) 2007-2009 Matthew West
# Licensed under the GNU General Public License version 2 or (at your
# option) any later version. See the file COPYING for details.

import os, sys, pyx
from numpy import *
sys.path.append("../tool")
from pmc_pyx import *


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

data_prefix = "aging_data/8"
coag_suffix = "nc"

smooth_window_len = 60

level = 4

grey_level = 0.8

aged_array = loadtxt("%s/aging_aged_%s" % (data_prefix, coag_suffix), unpack = True)
fresh_array = loadtxt("%s/aging_fresh_%s" % (data_prefix, coag_suffix), unpack = True)
emission_aged_array = loadtxt("%s/aging_emission_aged_%s" % (data_prefix, coag_suffix), unpack = True)
emission_fresh_array = loadtxt("%s/aging_emission_fresh_%s" % (data_prefix, coag_suffix), unpack = True)
loss_aged_array = loadtxt("%s/aging_loss_aged_%s" % (data_prefix, coag_suffix), unpack = True)
loss_fresh_array = loadtxt("%s/aging_loss_fresh_%s" % (data_prefix, coag_suffix), unpack = True)
transfer_to_aged_array = loadtxt("%s/aging_transfer_to_aged_%s" % (data_prefix, coag_suffix), unpack = True)
transfer_to_fresh_array = loadtxt("%s/aging_transfer_to_fresh_%s" % (data_prefix, coag_suffix), unpack = True)
height_array = loadtxt("%s/aging_height_%s" % (data_prefix, coag_suffix), unpack = True)

#outf = open("aging_k_nc", "w")

time = aged_array[0,:]
max_time_min = max(time) / 60
start_time_of_day_min = 6 * 60
height = height_array[1,:]
#time.tofile(outf, sep = " ", format = "%20f")
#outf.write("\n")
#for level in range(1, aged_array.shape[0]):

aged = aged_array[level,:]
fresh = fresh_array[level,:]

emission_aged = emission_aged_array[level,:]
emission_fresh = emission_fresh_array[level,:]

loss_aged = loss_aged_array[level,:]
loss_fresh = loss_fresh_array[level,:]

transfer_to_aged = transfer_to_aged_array[level,:]
transfer_to_fresh = transfer_to_fresh_array[level,:]

# hack due to halving
discard_threshold = 10
print "discard_threshold =", discard_threshold
median_loss_aged = median(loss_aged)
print "max(loss_aged) / median(loss_aged) =", max(loss_aged) / median_loss_aged
median_loss_fresh = median(loss_fresh)
print "max(loss_fresh) / median(loss_fresh) =", max(loss_fresh) / median_loss_fresh
for i in range(len(loss_aged)):
    if loss_aged[i] > discard_threshold * median_loss_aged:
        print "discarding loss_aged[%d]" % i
        loss_aged[i] = loss_aged[i - 1]
    if loss_fresh[i] > discard_threshold * median_loss_fresh:
        print "discarding loss_fresh[%d]" % i
        loss_fresh[i] = loss_fresh[i - 1]

transfer_net = transfer_to_aged - transfer_to_fresh

aged_raw = aged
fresh_raw = fresh
aged = smooth(aged, window_len = smooth_window_len)
fresh = smooth(fresh, window_len = smooth_window_len)

height_raw = height
height = smooth(height, window_len = smooth_window_len)

def delta(arr):
    return (arr[1:] - arr[:-1])

def mid(arr):
    return (arr[1:] + arr[:-1]) / 2.0

total = aged + fresh
aged_dot = delta(aged) / delta(time)
fresh_dot = delta(fresh) / delta(time)
height_dot = delta(height) / delta(time)
time_mid = mid(time)
aged_mid = mid(aged)
fresh_mid = mid(fresh)
height_mid = mid(height)

aged_raw_mid = mid(aged_raw)
fresh_raw_mid = mid(fresh_raw)
height_raw_mid = mid(height_raw)

aged_frac = aged_mid / (aged_mid + fresh_mid)
aged_frac_raw = aged_raw_mid / (aged_raw_mid + fresh_raw_mid)

emission_aged_end = emission_aged[1:]
emission_aged_end_raw = emission_aged_end
emission_aged_end = smooth(emission_aged_end, window_len = smooth_window_len)
emission_aged_end_rate = emission_aged_end / delta(time)
emission_aged_end_raw_rate = emission_aged_end_raw / delta(time)
emission_fresh_end = emission_fresh[1:]
emission_fresh_end_raw = emission_fresh_end
emission_fresh_end = smooth(emission_fresh_end, window_len = smooth_window_len)
emission_fresh_end_rate = emission_fresh_end / delta(time)
emission_fresh_end_raw_rate = emission_fresh_end_raw / delta(time)

loss_aged_end = loss_aged[1:]
loss_aged_end_raw = loss_aged_end
loss_aged_end = smooth(loss_aged_end, window_len = smooth_window_len)
loss_aged_end_rate = loss_aged_end / delta(time)
loss_aged_end_raw_rate = loss_aged_end_raw / delta(time)
loss_fresh_end = loss_fresh[1:]
loss_fresh_end_raw = loss_fresh_end
loss_fresh_end = smooth(loss_fresh_end, window_len = smooth_window_len)
loss_fresh_end_rate = loss_fresh_end / delta(time)
loss_fresh_end_raw_rate = loss_fresh_end_raw / delta(time)

transfer_to_aged_end = transfer_to_aged[1:]
transfer_to_aged_end_raw = transfer_to_aged_end
transfer_to_aged_end = smooth(transfer_to_aged_end, window_len = smooth_window_len)
transfer_to_aged_end_rate = transfer_to_aged_end / delta(time)
transfer_to_aged_end_raw_rate = transfer_to_aged_end_raw / delta(time)

transfer_to_fresh_end = transfer_to_fresh[1:]
transfer_to_fresh_end_raw = transfer_to_fresh_end
transfer_to_fresh_end = smooth(transfer_to_fresh_end, window_len = smooth_window_len)
transfer_to_fresh_end_rate = transfer_to_fresh_end / delta(time)
transfer_to_fresh_end_raw_rate = transfer_to_fresh_end_raw / delta(time)

transfer_net_end = transfer_net[1:]
transfer_net_end_raw = transfer_net_end
transfer_net_end = smooth(transfer_net_end, window_len = smooth_window_len)
transfer_net_end_rate = transfer_net_end / delta(time)
transfer_net_end_raw_rate = transfer_net_end_raw / delta(time)

dilution_rate_eff = dilution_rate + maximum(0, height_dot / height_mid)
loss_aged_dilution = dilution_rate_eff * aged_mid
loss_fresh_dilution = dilution_rate_eff * fresh_mid

#k_aged = (aged_dot + dilution_rate_eff * aged_mid - emission_aged_end_rate) / fresh_mid
#k_fresh = -(fresh_dot + dilution_rate_eff * fresh_mid - emission_fresh_end_rate) / fresh_mid
k_aged = (aged_dot - emission_aged_end_rate + loss_aged_end_rate) / fresh_mid
k_fresh = -(fresh_dot - emission_fresh_end_rate + loss_fresh_end_rate) / fresh_mid
k_transfer = transfer_net_end_rate / fresh_mid
k_transfer_raw = transfer_net_end_raw_rate / fresh_raw_mid

total_dot = delta(total) / delta(time)
total_mid = mid(total)
emission_total_end_rate = emission_aged_end_rate + emission_fresh_end_rate
loss_total_end_rate = loss_aged_end_rate + loss_fresh_end_rate
total_rhs = emission_total_end_rate - loss_total_end_rate
#total_rhs = -dilution_rate_eff * total_mid + emission_total_end_rate

tau_aged = 1.0 / k_aged
tau_fresh = 1.0 / k_fresh
tau_transfer = 1.0 / k_transfer
tau_transfer_raw = 1.0 / k_transfer_raw

######################################################################

g = pyx.graph.graphxy(
    width = 10,
    x = pyx.graph.axis.linear(title = r"time (s)",
                          painter = grid_painter),
    y = pyx.graph.axis.linear(painter = grid_painter),
    key = pyx.graph.key.key(pos = "tr"))

g.plot(
    pyx.graph.data.points(zip(time_mid, total_dot), x = 1, y = 2,
                          title = "total dot"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.red])])

g.plot(
    pyx.graph.data.points(zip(time_mid, total_rhs), x = 1, y = 2,
                          title = "total rhs"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.blue])])

g.writePDFfile("%s/aging_total_%d_%s.pdf" % (data_prefix, level, coag_suffix))

######################################################################

g = pyx.graph.graphxy(
    width = 10,
    x = pyx.graph.axis.linear(title = r"time (s)",
                          painter = grid_painter),
    y = pyx.graph.axis.linear(min = -0.5,
                              max = 0.5, #5,
                              title = r"$k$ ($\rm hour^{-1}$)",
                              painter = grid_painter),
    key = pyx.graph.key.key(pos = "tr"))

g.plot(
    pyx.graph.data.points(zip(time_mid, k_aged * 3600), x = 1, y = 2,
                          title = "k aged"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.red])])

g.plot(
    pyx.graph.data.points(zip(time_mid, k_fresh * 3600), x = 1, y = 2,
                          title = "k fresh"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.blue])])

g.plot(
    pyx.graph.data.points(zip(time_mid, k_transfer * 3600), x = 1, y = 2,
                          title = "k transfer"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.black])])

g.writePDFfile("%s/aging_k_%d_%s.pdf" % (data_prefix, level, coag_suffix))

######################################################################

g = pyx.graph.graphxy(
    width = 10,
    x = pyx.graph.axis.linear(title = r"time (s)",
                          painter = grid_painter),
    y = pyx.graph.axis.linear(min = -100,
                              max = 100,
                              painter = grid_painter),
    key = pyx.graph.key.key(pos = "tr"))

g.plot(
    pyx.graph.data.points(zip(time_mid, tau_aged / 3600), x = 1, y = 2,
                          title = "tau aged"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.red])])

g.plot(
    pyx.graph.data.points(zip(time_mid, tau_fresh / 3600), x = 1, y = 2,
                          title = "tau fresh"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.blue])])

g.plot(
    pyx.graph.data.points(zip(time_mid, tau_transfer / 3600), x = 1, y = 2,
                          title = "tau transfer"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.black])])

g.writePDFfile("%s/aging_tau_%d_%s.pdf" % (data_prefix, level, coag_suffix))

######################################################################

g = pyx.graph.graphxy(
    width = 10,
    x = pyx.graph.axis.linear(title = r"time (s)",
                          painter = grid_painter),
    y = pyx.graph.axis.linear(painter = grid_painter),
    key = pyx.graph.key.key(pos = "tr"))

g.plot(
    pyx.graph.data.points(zip(time_mid, emission_aged_end_raw_rate), x = 1, y = 2,
                          title = "emission aged raw"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.black])])

g.plot(
    pyx.graph.data.points(zip(time_mid, emission_fresh_end_raw_rate), x = 1, y = 2,
                          title = "emission fresh raw"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.green])])

g.plot(
    pyx.graph.data.points(zip(time_mid, emission_aged_end_rate), x = 1, y = 2,
                          title = "emission aged"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.red])])

g.plot(
    pyx.graph.data.points(zip(time_mid, emission_fresh_end_rate), x = 1, y = 2,
                          title = "emission fresh"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.blue])])

g.writePDFfile("%s/aging_emissions_%d_%s.pdf" % (data_prefix, level, coag_suffix))

######################################################################

g = pyx.graph.graphxy(
    width = 10,
    x = pyx.graph.axis.linear(title = r"time (s)",
                          painter = grid_painter),
    y = pyx.graph.axis.linear(painter = grid_painter),
    key = pyx.graph.key.key(pos = "tr"))

g.plot(
    pyx.graph.data.points(zip(time_mid, loss_aged_end_raw_rate), x = 1, y = 2,
                          title = "loss aged raw"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.black])])

g.plot(
    pyx.graph.data.points(zip(time_mid, loss_fresh_end_raw_rate), x = 1, y = 2,
                          title = "loss fresh raw"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.green])])

g.plot(
    pyx.graph.data.points(zip(time_mid, loss_aged_end_rate), x = 1, y = 2,
                          title = "loss aged"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.red])])

g.plot(
    pyx.graph.data.points(zip(time_mid, loss_fresh_end_rate), x = 1, y = 2,
                          title = "loss fresh"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.blue])])

g.plot(
    pyx.graph.data.points(zip(time_mid, loss_aged_dilution), x = 1, y = 2,
                          title = "dilution loss aged"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb(1,1,0)])])

g.plot(
    pyx.graph.data.points(zip(time_mid, loss_fresh_dilution), x = 1, y = 2,
                          title = "dilution loss fresh"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb(1,0,1)])])

g.writePDFfile("%s/aging_loss_%d_%s.pdf" % (data_prefix, level, coag_suffix))

######################################################################

g = pyx.graph.graphxy(
    width = 10,
    x = pyx.graph.axis.linear(title = r"time (s)",
                          painter = grid_painter),
    y = pyx.graph.axis.linear(painter = grid_painter),
    key = pyx.graph.key.key(pos = "tr"))

g.plot(
    pyx.graph.data.points(zip(time_mid, transfer_to_aged_end_raw_rate), x = 1, y = 2,
                          title = "transfer to aged raw"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.black])])

g.plot(
    pyx.graph.data.points(zip(time_mid, transfer_to_aged_end_rate), x = 1, y = 2,
                          title = "transfer to aged"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.red])])

g.plot(
    pyx.graph.data.points(zip(time_mid, transfer_to_fresh_end_raw_rate), x = 1, y = 2,
                          title = "transfer to fresh raw"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.green])])

g.plot(
    pyx.graph.data.points(zip(time_mid, transfer_to_fresh_end_rate), x = 1, y = 2,
                          title = "transfer to fresh"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.blue])])

g.writePDFfile("%s/aging_transfer_%d_%s.pdf" % (data_prefix, level, coag_suffix))

######################################################################

g = pyx.graph.graphxy(
    width = 10,
    x = pyx.graph.axis.linear(title = r"time (s)",
                          painter = grid_painter),
    y = pyx.graph.axis.linear(painter = grid_painter),
    key = pyx.graph.key.key(pos = "tr"))

g.plot(
    pyx.graph.data.points(zip(time_mid, aged_raw_mid), x = 1, y = 2,
                          title = "aged raw"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.black])])

g.plot(
    pyx.graph.data.points(zip(time_mid, fresh_raw_mid), x = 1, y = 2,
                          title = "fresh raw"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.green])])

g.plot(
    pyx.graph.data.points(zip(time_mid, aged_mid), x = 1, y = 2,
                          title = "aged"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.red])])

g.plot(
    pyx.graph.data.points(zip(time_mid, fresh_mid), x = 1, y = 2,
                          title = "fresh"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.blue])])

g.writePDFfile("%s/aging_raw_%d_%s.pdf" % (data_prefix, level, coag_suffix))

######################################################################

g = pyx.graph.graphxy(
    width = 10,
    x = graph.axis.linear(min = 0.,
                          max = max_time_min,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time
                                               = start_time_of_day_min),
                          title = "local standard time (LST) (hours:minutes)",
                          painter = grid_painter),
    y = pyx.graph.axis.linear(min = 0,
                              max = 100,
                              title = r"aged fraction (\%)",
                              painter = grid_painter))

g.plot(
    pyx.graph.data.points(zip(time_mid / 60, aged_frac_raw * 100), x = 1, y = 2),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.grey(grey_level)])])

g.plot(
    pyx.graph.data.points(zip(time_mid / 60, aged_frac * 100), x = 1, y = 2),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.black])])

g.writePDFfile("%s/aging_aged_frac_%d_%s.pdf" % (data_prefix, level, coag_suffix))

######################################################################

c = pyx.canvas.canvas()

g2 = c.insert(pyx.graph.graphxy(
    width = 10,
    x = graph.axis.linear(min = 0.,
                          max = max_time_min,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time
                                               = start_time_of_day_min),
                          title = "local standard time (LST) (hours:minutes)",
                          painter = grid_painter),
#    x = pyx.graph.axis.linear(title = r"time (s)",
#                          painter = grid_painter),
    y = pyx.graph.axis.log(reverse = 1,
                           min = 1e-1,
                           max = 1e5,
                           texter = pyx.graph.axis.texter.decimal(prefix = "-"),
                           painter = grid_painter)))

g1 = c.insert(pyx.graph.graphxy(
    width = 10,
    ypos = g2.height + 0.5,
    x = graph.axis.linkedaxis(g2.axes["x"],
                              painter = linked_grid_painter),
    y = pyx.graph.axis.log(min = 1e-1,
                           max = 1e5,
                              painter = grid_painter)))

gs = pyx.graph.graphxy(
    width = 10,
    x = graph.axis.linear(min = 0.,
                          max = max_time_min,
                          parter = graph.axis.parter.linear(tickdists
                                                            = [6 * 60, 3 * 60]),
                          texter = time_of_day(base_time
                                               = start_time_of_day_min),
                          title = "local standard time (LST) (hours:minutes)",
                          painter = grid_painter),
#    x = pyx.graph.axis.linear(title = r"time (s)",
#                              painter = grid_painter),
    y = pyx.graph.axis.log(min = 1e-1,
                           max = 1e5,
                           title = r"aging timescale $\tau$ (hours)",
                           painter = grid_painter))

def sign(x):
    if x > 0:
        return 1
    elif x < 0:
        return -1
    return 0

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

plot_data = zip(time_mid / 60, tau_transfer_raw / 3600)
chopped_data = chop_sign_data(plot_data)
for signed_data in chopped_data:
    if signed_data[0][1] > 0:
        g1.plot(
            pyx.graph.data.points(signed_data, x = 1, y = 2),
            styles = [pyx.graph.style.line(lineattrs = [pyx.color.grey(grey_level)])])
        gs.plot(
            pyx.graph.data.points(signed_data, x = 1, y = 2),
            styles = [pyx.graph.style.line(lineattrs = [pyx.color.grey(grey_level)])])
    else:
        signed_data = [[t,-d] for [t,d] in signed_data]
        g2.plot(
            pyx.graph.data.points(signed_data, x = 1, y = 2),
            styles = [pyx.graph.style.line(lineattrs = [pyx.color.grey(grey_level)])])
        gs.plot(
            pyx.graph.data.points(signed_data, x = 1, y = 2),
            styles = [pyx.graph.style.line(lineattrs = [pyx.color.grey(grey_level),
                                                        pyx.style.linestyle.dashed])])

plot_data = zip(time_mid / 60, tau_aged / 3600)
chopped_data = chop_sign_data(plot_data)
for signed_data in chopped_data:
    if signed_data[0][1] > 0:
        g1.plot(
            pyx.graph.data.points(signed_data, x = 1, y = 2),
            styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.red])])
#        gs.plot(
#            pyx.graph.data.points(signed_data, x = 1, y = 2),
#            styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.red])])
    else:
        signed_data = [[t,-d] for [t,d] in signed_data]
        g2.plot(
            pyx.graph.data.points(signed_data, x = 1, y = 2),
            styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.red])])
#        gs.plot(
#            pyx.graph.data.points(signed_data, x = 1, y = 2),
#            styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.red,
#                                                        pyx.style.linestyle.dashed])])

plot_data = zip(time_mid / 60, tau_fresh / 3600)
chopped_data = chop_sign_data(plot_data)
for signed_data in chopped_data:
    if signed_data[0][1] > 0:
        g1.plot(
            pyx.graph.data.points(signed_data, x = 1, y = 2),
            styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.blue])])
#        gs.plot(
#            pyx.graph.data.points(signed_data, x = 1, y = 2),
#            styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.blue])])
    else:
        signed_data = [[t,-d] for [t,d] in signed_data]
        g2.plot(
            pyx.graph.data.points(signed_data, x = 1, y = 2),
            styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.blue])])
#        gs.plot(
#            pyx.graph.data.points(signed_data, x = 1, y = 2),
#            styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.blue,
#                                                        pyx.style.linestyle.dashed])])

plot_data = zip(time_mid / 60, tau_transfer / 3600)
chopped_data = chop_sign_data(plot_data)
for signed_data in chopped_data:
    if signed_data[0][1] > 0:
        g1.plot(
            pyx.graph.data.points(signed_data, x = 1, y = 2),
            styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.black])])
        gs.plot(
            pyx.graph.data.points(signed_data, x = 1, y = 2),
            styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.black])])
    else:
        signed_data = [[t,-d] for [t,d] in signed_data]
        g2.plot(
            pyx.graph.data.points(signed_data, x = 1, y = 2),
            styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.black])])
        gs.plot(
            pyx.graph.data.points(signed_data, x = 1, y = 2),
            styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.black,
                                                        pyx.style.linestyle.dashed])])

c.writePDFfile("%s/aging_tau_log_%d_%s.pdf" % (data_prefix, level, coag_suffix))
gs.writePDFfile("%s/aging_tau_log_single_%d_%s.pdf" % (data_prefix, level, coag_suffix))

######################################################################
