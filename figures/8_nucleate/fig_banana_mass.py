#!/usr/bin/env python

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc

size_dist_array = np.loadtxt("data/banana_mass_dist_wc.txt")
times = np.loadtxt("data/banana_times_wc.txt")
diam = np.loadtxt("data/banana_diam_wc.txt")

mask = np.ma.make_mask(size_dist_array <= 0.0)
size_dist_array = np.ma.array(size_dist_array, mask = mask)
times = times / 3600.

plt.clf()
plt.pcolor(times, diam, size_dist_array.transpose(),norm = matplotlib.colors.LogNorm(), linewidths = 0.1)
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
plt.axis([times.min(), times.max(), diam.min(), diam.max()])
plt.xlabel("time (h)")
plt.ylabel("dry diameter (m)")
cbar = plt.colorbar()
cbar.set_label("number density (m^{-3})")
fig = plt.gcf()
fig.savefig("figs/banana_mass_wc.pdf")




