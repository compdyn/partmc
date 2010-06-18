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

num_avg_overall = np.zeros([12,config.i_loop_max])
mass_avg_overall = np.zeros([12,config.i_loop_max])
num_std_overall = np.zeros([12,config.i_loop_max])
mass_std_overall = np.zeros([12,config.i_loop_max])
i_counter = 0
for counter in ["10K_wei\\+1", "10K_flat", "10K_wei-1", "10K_wei-2", "10K_wei-3", "10K_wei-4","1K_wei\\+1", "1K_flat", "1K_wei-1", "1K_wei-2", "1K_wei-3", "1K_wei-4"]:
    f1 = "data/ensemble_number_%s.txt" % counter
    f2 = "data/ensemble_mass_%s.txt" % counter
    f3 = "data/ensemble_number_std_%s.txt" % counter
    f4 = "data/ensemble_mass_std_%s.txt" % counter

    num_avg_overall[i_counter,:] = np.loadtxt(f1)
    mass_avg_overall[i_counter,:] = np.loadtxt(f2)
    num_std_overall[i_counter,:] = np.loadtxt(f3)
    mass_std_overall[i_counter,:] = np.loadtxt(f4)

    i_counter +=1

x_array = [1, 10]

print 'check ', x_array[0],x_array[1]

plt.clf()
plt.errorbar(x_array[0], num_avg_overall[6,99], num_std_overall[6,99], marker='s', mfc='b', fmt='o', ecolor='g')
plt.errorbar(x_array[1], num_avg_overall[0,99], num_std_overall[0,99], marker='s', mfc='b',fmt='o', ecolor='g')
a = plt.gca()
a.set_xscale("log")
a.set_yscale("linear")
plt.axis([0.1, 100, 0.95e10, 1.05e10])
plt.grid(True)
plt.title("Number, alpha = 1")
plt.xlabel("particle number ")
plt.ylabel("mean number in m^{-3}")
fig = plt.gcf()
fig.savefig("figs/num_mean_v_particlen_wei+1.pdf")

plt.clf()
plt.errorbar(x_array[0], num_avg_overall[7,99], num_std_overall[7,99], marker='s', mfc='b',fmt='o', ecolor='g')
plt.errorbar(x_array[1], num_avg_overall[1,99], num_std_overall[1,99], marker='s', mfc='b',fmt='o', ecolor='g')
a = plt.gca()
a.set_xscale("log")
a.set_yscale("linear")
plt.axis([0.1, 100, 0.95e10, 1.05e10])
plt.grid(True)
plt.title("Number, alpha = 0")
plt.xlabel("particle number ")
plt.ylabel("mean number in m^{-3}")
fig = plt.gcf()
fig.savefig("figs/num_mean_v_particlen_flat.pdf")

plt.clf()
plt.errorbar(x_array[0], num_avg_overall[8,99], num_std_overall[8,99], marker='s', mfc='b',fmt='o', ecolor='g')
plt.errorbar(x_array[1], num_avg_overall[2,99], num_std_overall[2,99], marker='s', mfc='b',fmt='o', ecolor='g')
a = plt.gca()
a.set_xscale("log")
a.set_yscale("linear")
plt.axis([0.1, 100, 0.95e10, 1.05e10])
plt.grid(True)
plt.title("Number, alpha = -1")
plt.xlabel("particle number ")
plt.ylabel("mean number in m^{-3}")
fig = plt.gcf()
fig.savefig("figs/num_mean_v_particlen_wei-1.pdf")

plt.clf()
plt.errorbar(x_array[0], num_avg_overall[9,99], num_std_overall[9,99], marker='s', mfc='b',fmt='o', ecolor='g')
plt.errorbar(x_array[1], num_avg_overall[3,99], num_std_overall[3,99], marker='s', mfc='b',fmt='o', ecolor='g')
a = plt.gca()
a.set_xscale("log")
a.set_yscale("linear")
plt.axis([0.1, 100, 0.9e10, 1.1e10])
plt.grid(True)
plt.title("Number, alpha = -2")
plt.xlabel("particle number ")
plt.ylabel("mean number in m^{-3}")
fig = plt.gcf()
fig.savefig("figs/num_mean_v_particlen_wei-2.pdf")

plt.clf()
plt.errorbar(x_array[0], num_avg_overall[10,99], num_std_overall[10,99], marker='s', mfc='b',fmt='o', ecolor='g')
plt.errorbar(x_array[1], num_avg_overall[4,99], num_std_overall[4,99], marker='s', mfc='b',fmt='o', ecolor='g')
a = plt.gca()
a.set_xscale("log")
a.set_yscale("linear")
plt.axis([0.1, 100, 0.5e10, 1.5e10])
plt.grid(True)
plt.title("Number, alpha = -3")
plt.xlabel("particle number ")
plt.ylabel("mean number in m^{-3}")
fig = plt.gcf()
fig.savefig("figs/num_mean_v_particlen_wei-3.pdf")

plt.clf()
plt.errorbar(x_array[0], num_avg_overall[11,99], num_std_overall[11,99], marker='s', mfc='b',fmt='o', ecolor='g')
plt.errorbar(x_array[1], num_avg_overall[5,99], num_std_overall[5,99], marker='s', mfc='b',fmt='o', ecolor='g')
a = plt.gca()
a.set_xscale("log")
a.set_yscale("linear")
plt.axis([0.1, 100, 0.4e10, 1.5e10])
plt.grid(True)
plt.title("Number, alpha = -4")
plt.xlabel("particle number ")
plt.ylabel("mean number in m^{-3}")
fig = plt.gcf()
fig.savefig("figs/num_mean_v_particlen_wei-4.pdf")


plt.clf()
plt.errorbar(x_array[0], mass_avg_overall[6,99], mass_std_overall[6,99], marker='s', mfc='b', fmt='o', ecolor='g')
plt.errorbar(x_array[1], mass_avg_overall[0,99], mass_std_overall[0,99], marker='s', mfc='b',fmt='o', ecolor='g')
a = plt.gca()
a.set_xscale("log")
a.set_yscale("linear")
plt.axis([0.1, 100, 20.5, 21.5])
plt.xlabel("particle mass ")
plt.ylabel("mean mass in \mu m^{-3}")
plt.title("Mass, alpha = 1")
fig = plt.gcf()
fig.savefig("figs/mass_mean_v_particlen_wei+1.pdf")

plt.clf()
plt.errorbar(x_array[0], mass_avg_overall[7,99], mass_std_overall[7,99], marker='s', mfc='b',fmt='o', ecolor='g')
plt.errorbar(x_array[1], mass_avg_overall[1,99], mass_std_overall[1,99], marker='s', mfc='b',fmt='o', ecolor='g')
a = plt.gca()
a.set_xscale("log")
a.set_yscale("linear")
plt.axis([0.1, 100, 20.5, 21.5])
plt.xlabel("particle mass ")
plt.ylabel("mean mass in \mu m^{-3}")
plt.title("Mass, alpha = 0")
fig = plt.gcf()
fig.savefig("figs/mass_mean_v_particlen_flat.pdf")

plt.clf()
plt.errorbar(x_array[0], mass_avg_overall[8,99], mass_std_overall[8,99], marker='s', mfc='b',fmt='o', ecolor='g')
plt.errorbar(x_array[1], mass_avg_overall[2,99], mass_std_overall[2,99], marker='s', mfc='b',fmt='o', ecolor='g')
a = plt.gca()
a.set_xscale("log")
a.set_yscale("linear")
plt.axis([0.1, 100, 20.5, 21.5])
plt.xlabel("particle mass ")
plt.ylabel("mean mass in \mu m^{-3}")
plt.title("Mass, alpha = -1")
fig = plt.gcf()
fig.savefig("figs/mass_mean_v_particlen_wei-1.pdf")

plt.clf()
plt.errorbar(x_array[0], mass_avg_overall[9,99], mass_std_overall[9,99], marker='s', mfc='b',fmt='o', ecolor='g')
plt.errorbar(x_array[1], mass_avg_overall[3,99], mass_std_overall[3,99], marker='s', mfc='b',fmt='o', ecolor='g')
a = plt.gca()
a.set_xscale("log")
a.set_yscale("linear")
plt.axis([0.1, 100, 20.5, 21.5])
plt.xlabel("particle mass ")
plt.ylabel("mean mass in \mu m^{-3}")
plt.title("Mass, alpha = -2")
fig = plt.gcf()
fig.savefig("figs/mass_mean_v_particlen_wei-2.pdf")

plt.clf()
plt.errorbar(x_array[0], mass_avg_overall[10,99], mass_std_overall[10,99], marker='s', mfc='b',fmt='o', ecolor='g')
plt.errorbar(x_array[1], mass_avg_overall[4,99], mass_std_overall[4,99], marker='s', mfc='b',fmt='o', ecolor='g')
a = plt.gca()
a.set_xscale("log")
a.set_yscale("linear")
plt.axis([0.1, 100, 20.5, 21.5])
plt.xlabel("particle mass ")
plt.ylabel("mean mass in \mu m^{-3}")
plt.title("Mass, alpha = -3")
fig = plt.gcf()
fig.savefig("figs/mass_mean_v_particlen_wei-3.pdf")


plt.clf()
plt.errorbar(x_array[0], mass_avg_overall[11,99], mass_std_overall[11,99], marker='s', mfc='b',fmt='o', ecolor='g')
plt.errorbar(x_array[1], mass_avg_overall[5,99], mass_std_overall[5,99], marker='s', mfc='b',fmt='o', ecolor='g')
a = plt.gca()
a.set_xscale("log")
a.set_yscale("linear")
plt.axis([0.1, 100, 20.5, 21.5])
plt.xlabel("particle mass ")
plt.ylabel("mean mass in \mu m^{-3}")
plt.title("Mass, alpha = -4")
fig = plt.gcf()
fig.savefig("figs/mass_mean_v_particlen_wei-4.pdf")







