#!/usr/bin/env python

import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt

data = np.loadtxt("urban_plume_32_none_dist/diam_bc_num_0001_00000008_error.txt")
plt.loglog(data[:,0], data[:,1], '-D', label="Prob mix = 0")
plt.hold(True)
data = np.loadtxt("urban_plume_32_local1_dist/diam_bc_num_0001_00000008_error.txt")
plt.loglog(data[:,0], data[:,1], '-o', label="Prob mix = 0.01")
data = np.loadtxt("urban_plume_32_local2_dist/diam_bc_num_0001_00000008_error.txt")
plt.loglog(data[:,0], data[:,1], '-x', label="Prob mix = 0.1")
data = np.loadtxt("urban_plume_32_local3_dist/diam_bc_num_0001_00000008_error.txt")
plt.loglog(data[:,0], data[:,1], '-s', label="Prob mix = 1")

x0 = 1.0
y0 = 2.0e5
x1 = x0 / 4
y1 = y0 * 2
x2 = x0 * 10000
y2 = y0 / 100
plt.loglog([x1, x2], [y1, y2], '-', label="slope -0.5")

plt.xlabel("cores")
plt.ylabel("RMSE for D/BC number conc. at t = 7 hr")
plt.title("32 particles per core, 1441 timesteps, 60 s timestep")
plt.xlim([1, 1e3])
plt.ylim([8e3, 3e5])
plt.grid(True)
plt.legend(loc="upper right")
fig = plt.gcf()
fig.savefig("parallel_mixing_32_error.pdf")
