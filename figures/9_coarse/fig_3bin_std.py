#!/usr/bin/env python

import scipy.io
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc
import config

array_error_num_10K = np.loadtxt("data/num_1d_3bin_10K.txt")
array_error_num_100K = np.loadtxt("data/num_1d_3bin_100K.txt")
array_error_mass_10K = np.loadtxt("data/mass_1d_3bin_10K.txt")
array_error_mass_100K = np.loadtxt("data/mass_1d_3bin_100K.txt")

x_array = [0, -1, -2, -3, -4]

plt.clf()
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
plt.plot(x_array, array_error_num_100K[:,0], 'r-', label = 'number 100K')
plt.plot(x_array, array_error_num_10K[:,0], 'b-', label = 'number 10K')
plt.plot(x_array, array_error_mass_100K[:,0], 'r--', label = 'mass 100K')
plt.plot(x_array, array_error_mass_10K[:,0], 'b--', label = 'mass 10K')
plt.legend(loc = 'upper right')
plt.xlabel("weighting exponent ")
plt.ylabel("averaged std/avg (small particles)")
fig = plt.gcf()
fig.savefig("figs/std_total_averages_small.pdf")

plt.clf()
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
plt.plot(x_array, array_error_num_100K[:,1], 'r-', label = 'number 100K')
plt.plot(x_array, array_error_num_10K[:,1], 'b-', label = 'number 10K')
plt.plot(x_array, array_error_mass_100K[:,1], 'r--', label = 'mass 100K')
plt.plot(x_array, array_error_mass_10K[:,1], 'b--', label = 'mass 10K')
plt.legend(loc = 'upper right')
plt.xlabel("weighting exponent ")
plt.ylabel("averaged std/avg (medium particles)")
fig = plt.gcf()
fig.savefig("figs/std_total_averages_med.pdf")

plt.clf()
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
plt.plot(x_array, array_error_num_100K[:,2], 'r-', label = 'number 100K')
plt.plot(x_array, array_error_num_10K[:,2], 'b-', label = 'number 10K')
plt.plot(x_array, array_error_mass_100K[:,2], 'r--', label = 'mass 100K')
plt.plot(x_array, array_error_mass_10K[:,2], 'b--', label = 'mass 10K')
plt.legend(loc = 'upper right')
plt.xlabel("weighting exponent ")
plt.ylabel("averaged std/avg (large particles)")
fig = plt.gcf()
fig.savefig("figs/std_total_averages_large.pdf")




    



