#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

a = np.loadtxt("data/maximum_ss.txt")
#b = np.loadtxt("data/ccn_cn_ratio_wc.txt")

time = range(0,49)

print time,a[0,:]

plt.clf()
plt.plot(time,a[0,:], "*", label = "ref")
plt.plot(time,a[1,:], "*", label = "comp")
plt.plot(time,a[2,:], "*", label = "size")
plt.plot(time,a[3,:], "*", label = "both")
plt.legend(loc = 'lower right')
plt.grid(True)
plt.axis([0, 48, 0.1, 0.35])
fig = plt.gcf()
plt.xlabel("time (hour)")
plt.ylabel("maximum supersaturation (%)")
fig.savefig("figs/maximum_ss_wc.pdf")

#plt.clf()
#plt.plot(time,b[0,:], label = "ref")
#plt.plot(time,b[1,:], label = "comp")
#plt.plot(time,b[2,:], label = "size")
#plt.plot(time,b[3,:], label = "both")
#plt.legend(loc = 'upper left')
#fig = plt.gcf()
#plt.xlabel("time (hour)")
#plt.ylabel("ccn/cn")
#fig.savefig("figs/ccn_cn_ratio_wc.pdf")




