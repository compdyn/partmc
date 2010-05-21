#!/usr/bin/env python

import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt

data = np.loadtxt("times_32_none_dist.txt")
cores = data[:-3,0]
plt.loglog(cores, cores * data[0,1] / data[:-3,1], '-D', label="Prob mix = 0")
plt.hold(True)

data = np.loadtxt("times_local1_dist.txt")
cores = data[:,0] * 8
cores[0] = 1
plt.loglog(cores, cores * data[0,1] / data[:,1], '-o', label="Prob mix = 0.01")
plt.hold(True)
data = np.loadtxt("times_local2_dist.txt")
plt.loglog(cores, cores * data[0,1] / data[:,1], '-x', label="Prob mix = 0.1")
data = np.loadtxt("times_local3_dist.txt")
plt.loglog(cores, cores * data[0,1] / data[:,1], '-s', label="Prob mix = 1")

plt.loglog([1e-1, 1e4], [1e-1, 1e4], '-', label="slope 1")

plt.xlabel("cores")
plt.ylabel("speedup factor")
plt.title("1024 particles per core, 1441 timesteps, 60 s timestep")
plt.xlim([1, 1e3])
plt.ylim([1, 1e2])
plt.grid(True)
plt.legend(loc="lower right")
fig = plt.gcf()
fig.savefig("parallel_mixing_speedup.pdf")
