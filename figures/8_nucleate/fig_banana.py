#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
sys.path.append("../../tool")
import mpl_helper
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
import partmc

size_dist_array = np.loadtxt("data/banana_size_dist_wc_wei-1_lowbg2.txt")
times = np.loadtxt("data/banana_times_wc_wei-1_lowbg2.txt")
diam = np.loadtxt("data/banana_diam_wc_wei-1_lowbg2.txt")

mask = np.ma.make_mask(size_dist_array <= 0.0)
size_dist_array = np.ma.array(size_dist_array, mask = mask)
times = times / 3600.

(fig, axis) = mpl_helper.make_fig(figure_width = 4)

plt.pcolor(times, diam, size_dist_array.transpose(),norm = matplotlib.colors.LogNorm(vmin=1e1, vmax=1e6), linewidths = 0.1)
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
plt.axis([times.min(), times.max(), 1e-3, diam.max()])
plt.xlabel("time / h")
plt.ylabel(r"dry diameter / $\rm \mu m$")
plt.grid(True)
cbar = plt.colorbar(format=matplotlib.ticker.LogFormatterMathtext())
cbar.set_label(r"number concentration / $\rm cm^{-3}$")
fig = plt.gcf()
fig.savefig("figs/banana_wc_wei-1_lowbg2.pdf")




