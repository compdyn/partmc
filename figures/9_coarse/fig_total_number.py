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

time_array = np.loadtxt("data/time_array_10K_flat.txt")
array_num = np.loadtxt("data/array_number_10K_flat.txt")
num_avg = np.loadtxt("data/average_number_10K_flat.txt")
num_std = np.loadtxt("data/std_number_10K_flat.txt")
array_mass = np.loadtxt("data/array_mass_10K_flat.txt")
mass_avg = np.loadtxt("data/average_mass_10K_flat.txt")
mass_std = np.loadtxt("data/std_mass_10K_flat.txt")

plt.clf()
#for i_loop in range(0,config.i_loop_max):
#    plt.plot(time_array[:], array_num[:,i_loop], 'k')
plt.errorbar(time_array, num_avg, num_std)
plt.xlim(0, 24)
plt.xlabel("time ")
plt.ylabel("number concentration in m^{-3}")
plt.grid(True)
#plt.title("100K wei-3")
fig = plt.gcf()
fig.savefig("figs/number_10K_flat.pdf")

plt.clf()
#for i_loop in range(0,config.i_loop_max):
#    plt.plot(time_array, array_mass[:,i_loop], 'k')
plt.errorbar(time_array, mass_avg, mass_std)
plt.xlim(0, 24)
plt.xlabel("time ")
plt.ylabel("mass concentration in \mu g m^{-3}")
#plt.title("100K wei-3")
plt.grid(True)
fig = plt.gcf()
fig.savefig("figs/mass_10K_flat.pdf")

plt.clf()
plt.plot(time_array, mass_std/mass_avg)
plt.xlim(0, 24)
plt.xlabel("time ")
plt.ylabel("CV(M(\infty, t))")
#plt.title("100K wei-3")
plt.grid(True)
fig = plt.gcf()
fig.savefig("figs/std_mass_10K_flat.pdf")

plt.clf()
plt.plot(time_array, num_std/num_avg)
plt.xlim(0, 24)
plt.xlabel("time ")
plt.ylabel("CV(N(\infty, t))")
#plt.title("100K wei-3")
plt.grid(True)
fig = plt.gcf()
fig.savefig("figs/std_number_10K_flat.pdf")




    



