#!/usr/bin/env python

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc
import config

time_array = np.loadtxt("data/time_array.txt")
array_num = np.loadtxt("data/array_number.txt")
num_avg = np.loadtxt("data/average_number.txt")
num_std = np.loadtxt("data/std_number.txt")
array_mass = np.loadtxt("data/array_mass.txt")
mass_avg = np.loadtxt("data/average_mass.txt")
mass_std = np.loadtxt("data/std_mass.txt")

plt.clf()
for i_loop in range(0,config.i_loop_max):
    plt.plot(time_array[:], array_num[:,i_loop], 'k')
plt.errorbar(time_array, num_avg, num_std)
plt.xlabel("time ")
plt.ylabel("number concentration in m^{-3}")
plt.title("10K flat")
fig = plt.gcf()
fig.savefig("figs/number_10K_wei-3.pdf")

plt.clf()
for i_loop in range(0,config.i_loop_max):
    plt.plot(time_array, array_mass[:,i_loop], 'k')
plt.errorbar(time_array, mass_avg, mass_std)
plt.xlabel("time ")
plt.ylabel("mass concentration in \mu g m^{-3}")
plt.title("10K flat")
fig = plt.gcf()
fig.savefig("figs/mass_10K_wei-3.pdf")




    



