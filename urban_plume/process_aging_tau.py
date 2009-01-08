#!/usr/bin/env python
# Copyright (C) 2007, 2008 Matthew West
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

aged_array = loadtxt("aging_a_wc", unpack = True)
fresh_array = loadtxt("aging_f_wc", unpack = True)
emission_aged_array = loadtxt("aging_ea_wc", unpack = True)
emission_fresh_array = loadtxt("aging_ef_wc", unpack = True)
height_array = loadtxt("aging_h_wc", unpack = True)

#outf = open("aging_k_nc", "w")

time = aged_array[0,:]
height = height_array[1,:]
#time.tofile(outf, sep = " ", format = "%20f")
#outf.write("\n")
#for level in range(1, aged_array.shape[0]):

level = 10
aged = aged_array[level,:]
fresh = fresh_array[level,:]

emission_aged = emission_aged_array[level,:]
emission_fresh = emission_fresh_array[level,:]

#aged = smooth(aged)
#fresh = smooth(fresh)

total = aged + fresh
aged_dot = (aged[1:] - aged[:-1]) / (time[1:] - time[:-1])
fresh_dot = (fresh[1:] - fresh[:-1]) / (time[1:] - time[:-1])
height_dot = (height[1:] - height[:-1]) / (time[1:] - time[:-1])
time_mid = (time[1:] + time[:-1]) / 2.0
aged_mid = (aged[1:] + aged[:-1]) / 2.0
fresh_mid = (fresh[1:] + fresh[:-1]) / 2.0
height_mid = (height[1:] + height[:-1]) / 2.0
emission_aged_end = emission_aged[1:]
#emission_aged_end = smooth(emission_aged_end, window_len = 100)
emission_aged_end_rate = emission_aged_end / (time[1:] - time[:-1])
emission_fresh_end = emission_fresh[1:]
#emission_fresh_end = smooth(emission_fresh_end, window_len = 100)
emission_fresh_end_rate = emission_fresh_end / (time[1:] - time[:-1])

dilution_rate_eff = dilution_rate + maximum(0, height_dot / height_mid);

k_aged = (aged_dot + dilution_rate_eff * aged_mid - emission_aged_end_rate) / fresh_mid
k_fresh = -(fresh_dot + dilution_rate_eff * fresh_mid - emission_fresh_end_rate) / fresh_mid

total_dot = (total[1:] - total[:-1]) / (time[1:] - time[:-1])
total_mid = (total[1:] + total[:-1]) / 2.0
emission_total_end_rate = emission_aged_end_rate + emission_fresh_end_rate
total_rhs = -dilution_rate_eff * total_mid + emission_total_end_rate

#k_aged = smooth(k_aged)
#k_fresh = smooth(k_fresh)

tau_aged = 1.0 / k_aged
tau_fresh = 1.0 / k_fresh

#tau_aged = smooth(tau_aged, window_len = 100)
#tau_fresh = smooth(tau_fresh, window_len = 100)

#    k_aged.tofile(outf, sep = " ", format = "%20f")
#    outf.write("\n")
#    k_fresh.tofile(outf, sep = " ", format = "%20f")
#    outf.write("\n")

#outf.close()

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

g.writePDFfile("aging_total.pdf")

######################################################################

g = pyx.graph.graphxy(
    width = 10,
    x = pyx.graph.axis.linear(title = r"time (s)",
                          painter = grid_painter),
    y = pyx.graph.axis.linear(min = -0.5,
                              max = 5,
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

g.writePDFfile("aging_k_nc.pdf")

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

g.writePDFfile("aging_tau_nc.pdf")

######################################################################

g = pyx.graph.graphxy(
    width = 10,
    x = pyx.graph.axis.linear(title = r"time (s)",
                          painter = grid_painter),
    y = pyx.graph.axis.linear(painter = grid_painter),
    key = pyx.graph.key.key(pos = "tr"))

g.plot(
    pyx.graph.data.points(zip(time_mid, emission_aged_end_rate), x = 1, y = 2,
                          title = "emission aged"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.red])])

g.plot(
    pyx.graph.data.points(zip(time_mid, emission_fresh_end_rate), x = 1, y = 2,
                          title = "emission fresh"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.blue])])

g.writePDFfile("aging_emissions.pdf")

######################################################################

g = pyx.graph.graphxy(
    width = 10,
    x = pyx.graph.axis.linear(title = r"time (s)",
                          painter = grid_painter),
    y = pyx.graph.axis.linear(painter = grid_painter),
    key = pyx.graph.key.key(pos = "tr"))

g.plot(
    pyx.graph.data.points(zip(time_mid, aged_mid), x = 1, y = 2,
                          title = "aged"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.red])])

g.plot(
    pyx.graph.data.points(zip(time_mid, fresh_mid), x = 1, y = 2,
                          title = "fresh"),
    styles = [pyx.graph.style.line(lineattrs = [pyx.color.rgb.blue])])

g.writePDFfile("aging_raw.pdf")

######################################################################
