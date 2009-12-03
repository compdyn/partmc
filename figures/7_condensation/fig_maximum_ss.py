#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc

a = np.loadtxt("data/maximum_ss.txt")
b = np.loadtxt("data/ccn_cn_ratio.txt")

time = [1, 7, 15, 24]

plt.clf()
plt.plot(time,a[0,:], "*", label = "ref")
plt.plot(time,a[1,:], "*", label = "comp")
plt.plot(time,a[2,:], "*", label = "size")
plt.plot(time,a[3,:], "*", label = "both")
plt.legend(loc = 'upper right')
fig = plt.gcf()
plt.xlabel("time (hour)")
plt.ylabel("maximum supersaturation (%)")
fig.savefig("figs/maximum_ss.pdf")

plt.clf()
plt.plot(time,b[0,:], label = "ref")
plt.plot(time,b[1,:], label = "comp")
plt.plot(time,b[2,:], label = "size")
plt.plot(time,b[3,:], label = "both")
plt.legend(loc = 'upper left')
fig = plt.gcf()
plt.xlabel("time (hour)")
plt.ylabel("ccn/cn")
fig.savefig("figs/ccn_cn_ratio.pdf")




