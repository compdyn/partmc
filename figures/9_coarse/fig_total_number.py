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

time_array = np.loadtxt("data/time_array_100K_wei-3.txt")
array_num = np.loadtxt("data/array_number_100K_wei-3.txt")
num_avg = np.loadtxt("data/average_number_100K_wei-3.txt")
num_std = np.loadtxt("data/std_number_100K_wei-3.txt")
array_mass = np.loadtxt("data/array_mass_100K_wei-3.txt")
mass_avg = np.loadtxt("data/average_mass_100K_wei-3.txt")
mass_std = np.loadtxt("data/std_mass_100K_wei-3.txt")

plt.clf()
for i_loop in range(0,config.i_loop_max):
    plt.plot(time_array[:], array_num[:,i_loop], 'k')
plt.errorbar(time_array, num_avg, num_std)
plt.xlabel("time ")
plt.ylabel("number concentration in m^{-3}")
plt.title("100K wei-3")
fig = plt.gcf()
fig.savefig("figs/number_100K_wei-3.pdf")

plt.clf()
for i_loop in range(0,config.i_loop_max):
    plt.plot(time_array, array_mass[:,i_loop], 'k')
plt.errorbar(time_array, mass_avg, mass_std)
plt.xlabel("time ")
plt.ylabel("mass concentration in \mu g m^{-3}")
plt.title("100K wei-3")
fig = plt.gcf()
fig.savefig("figs/mass_100K_wei-3.pdf")

plt.clf()
plt.plot(time_array, mass_std/mass_avg)
plt.xlabel("time ")
plt.ylabel("rel. std in \mu g m^{-3}")
plt.title("100K wei-3")
fig = plt.gcf()
fig.savefig("figs/std_mass_100K_wei-3.pdf")

plt.clf()
plt.plot(time_array, num_std/num_avg)
plt.xlabel("time ")
plt.ylabel("rel. std in m^{-3}")
plt.title("100K wei-3")
fig = plt.gcf()
fig.savefig("figs/std_number_100K_wei-3.pdf")




    



