#!/usr/bin/env python

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

dist_array = np.loadtxt("data/banana_so4_dist.txt")
times = np.loadtxt("data/banana_times.txt")
diam = np.loadtxt("data/banana_diam.txt")

mask = np.ma.make_mask(dist_array <= 0.0)
dist_array = np.ma.array(dist_array, mask = mask)
times = times / 3600.

plt.clf()
plt.pcolor(times, diam, dist_array.transpose(),norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
plt.axis([times.min(), times.max(), diam.min(), diam.max()])
plt.xlabel("time (h)")
plt.ylabel("dry diameter (m)")
cbar = plt.colorbar()
cbar.set_label("mass density (kg m^{-3})")
fig = plt.gcf()
fig.savefig("figs/banana_so4_wc.pdf")




