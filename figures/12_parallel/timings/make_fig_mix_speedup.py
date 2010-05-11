#!/usr/bin/env python

import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt

data = np.loadtxt("times_local1_dist.txt")
cores = data[:,0] * 8
cores[0] = 1
plt.plot(cores, cores * data[0,1] / data[:,1], label="Prob mix = 0.01")
plt.hold(True)
data = np.loadtxt("times_local2_dist.txt")
plt.plot(cores, cores * data[0,1] / data[:,1], label="Prob mix = 0.1")
data = np.loadtxt("times_local3_dist.txt")
plt.plot(cores, cores * data[0,1] / data[:,1], label="Prob mix = 0.5")

plt.xlabel("cores")
plt.ylabel("speedup factor")
plt.title("1000 particles per core, 1441 timesteps")
plt.xlim([0, 128])
#plt.ylim([0, 900])
plt.grid(True)
plt.legend(loc="lower right")
fig = plt.gcf()
fig.savefig("parallel_mixing_speedup.pdf")
