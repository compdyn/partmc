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

num_avg_overall = np.zeros([14,config.i_loop_max])
mass_avg_overall = np.zeros([14,config.i_loop_max])
num_std_overall = np.zeros([14,config.i_loop_max])
mass_std_overall = np.zeros([14,config.i_loop_max])
i_counter = 0
for counter in ["10K_wei\\+1", "10K_flat", "10K_wei-1", "10K_wei-2", "10K_wei-3", "10K_wei-4","10K_mfa", "1K_wei\\+1", "1K_flat", "1K_wei-1", "1K_wei-2", "1K_wei-3", "1K_wei-4", "1K_mfa"]:
    f1 = "data/ensemble_number_%s.txt" % counter
    f2 = "data/ensemble_mass_%s.txt" % counter
    f3 = "data/ensemble_number_std_%s.txt" % counter
    f4 = "data/ensemble_mass_std_%s.txt" % counter

    num_avg_overall[i_counter,:] = np.loadtxt(f1)
    mass_avg_overall[i_counter,:] = np.loadtxt(f2)
    num_std_overall[i_counter,:] = np.loadtxt(f3)
    mass_std_overall[i_counter,:] = np.loadtxt(f4)

    i_counter +=1


x_array = np.arange(1,101)
print 'x_array ', x_array[1:], num_avg_overall[0,1:]

plt.clf()
plt.plot(x_array[1:], num_avg_overall[0,1:], 'c', label = '10K wei+1')
plt.plot(x_array[1:], num_avg_overall[1,1:], 'b', label = '10K flat')
plt.plot(x_array[1:], num_avg_overall[2,1:], 'g', label = '10K wei-1')
plt.plot(x_array[1:], num_avg_overall[3,1:], 'r', label = '10K wei-2')
plt.plot(x_array[1:], num_avg_overall[4,1:], 'k', label = '10K wei-3')
plt.plot(x_array[1:], num_avg_overall[5,1:], 'm', label = '10K wei-4')
plt.plot(x_array[1:], num_avg_overall[6,1:], 'y', label = '10K mfa')
plt.xlabel("ensemble size ")
plt.ylabel("average number concentration in m^{-3}")
plt.grid()
plt.legend(loc = 'lower right')
fig = plt.gcf()
fig.savefig("figs/number_ensem_10K.pdf")

plt.clf()
plt.plot(x_array[3:], num_std_overall[0,3:]/num_avg_overall[0,3:], 'c', label = '10K wei+1')
plt.plot(x_array[3:], num_std_overall[1,3:]/num_avg_overall[1,3:], 'b', label = '10K flat')
plt.plot(x_array[3:], num_std_overall[2,3:]/num_avg_overall[2,3:], 'g', label = '10K wei-1')
plt.plot(x_array[3:], num_std_overall[3,3:]/num_avg_overall[3,3:], 'r', label = '10K wei-2')
plt.plot(x_array[3:], num_std_overall[4,3:]/num_avg_overall[4,3:], 'k', label = '10K wei-3')
plt.plot(x_array[3:], num_std_overall[5,3:]/num_avg_overall[5,3:], 'm', label = '10K wei-4')
plt.plot(x_array[3:], num_std_overall[6,3:]/num_avg_overall[6,3:], 'y', label = '10K mfa')
plt.xlabel("ensemble size ")
plt.ylabel("average std number concentration")
plt.grid()
plt.legend(loc = 'lower right')
fig = plt.gcf()
fig.savefig("figs/number_std_ensem_10K.pdf")

plt.clf()
plt.plot(x_array[1:], mass_avg_overall[0,1:], 'c', label = '10K wei+1')
plt.plot(x_array[1:], mass_avg_overall[1,1:], 'b', label = '10K flat')
plt.plot(x_array[1:], mass_avg_overall[2,1:], 'g', label = '10K wei-1')
plt.plot(x_array[1:], mass_avg_overall[3,1:], 'r', label = '10K wei-2')
plt.plot(x_array[1:], mass_avg_overall[4,1:], 'k', label = '10K wei-3')
plt.plot(x_array[1:], mass_avg_overall[5,1:], 'm', label = '10K wei-4')
plt.plot(x_array[1:], mass_avg_overall[6,1:], 'y', label = '10K mfa')
plt.xlabel("ensemble size ")
plt.ylabel("average mass concentration in \mu g m^{-3}")
plt.grid()
plt.legend(loc = 'lower right')
fig = plt.gcf()
fig.savefig("figs/mass_ensem_10K.pdf")

plt.clf()
plt.plot(x_array[3:], mass_std_overall[0,3:]/mass_avg_overall[0,3:], 'c', label = '10K wei+1')
plt.plot(x_array[3:], mass_std_overall[1,3:]/mass_avg_overall[1,3:], 'b', label = '10K flat')
plt.plot(x_array[3:], mass_std_overall[2,3:]/mass_avg_overall[2,3:], 'g', label = '10K wei-1')
plt.plot(x_array[3:], mass_std_overall[3,3:]/mass_avg_overall[3,3:], 'r', label = '10K wei-2')
plt.plot(x_array[3:], mass_std_overall[4,3:]/mass_avg_overall[4,3:], 'k', label = '10K wei-3')
plt.plot(x_array[3:], mass_std_overall[5,3:]/mass_avg_overall[5,3:], 'm', label = '10K wei-4')
plt.plot(x_array[3:], mass_std_overall[6,3:]/mass_avg_overall[6,3:], 'y', label = '10K mfa')
plt.xlabel("ensemble size ")
plt.ylabel("average std mass concentration")
plt.grid()
plt.legend(loc = 'lower right')
fig = plt.gcf()
fig.savefig("figs/mass_std_ensem_10K.pdf")


plt.clf()
plt.plot(x_array[1:], num_avg_overall[7,1:], 'c', label = '1K wei+1')
plt.plot(x_array[1:], num_avg_overall[8,1:], 'b', label = '1K flat')
plt.plot(x_array[1:], num_avg_overall[9,1:], 'g', label = '1K wei-1')
plt.plot(x_array[1:], num_avg_overall[10,1:], 'r', label = '1K wei-2')
plt.plot(x_array[1:], num_avg_overall[11,1:], 'k', label = '1K wei-3')
plt.plot(x_array[1:], num_avg_overall[12,1:], 'm', label = '1K wei-4')
plt.plot(x_array[1:], num_avg_overall[13,1:], 'y', label = '1K mfa')
plt.xlabel("ensemble size ")
plt.ylabel("average number concentration in m^{-3}")
plt.grid()
plt.legend(loc = 'lower right')
fig = plt.gcf()
fig.savefig("figs/number_ensem_1K.pdf")
print 'check ', num_std_overall[6,3:],num_avg_overall[6,3:]
plt.clf()
plt.plot(x_array[3:], num_std_overall[7,3:]/num_avg_overall[7,3:], 'c', label = '1K wei+1')
plt.plot(x_array[3:], num_std_overall[8,3:]/num_avg_overall[8,3:], 'b', label = '1K flat')
plt.plot(x_array[3:], num_std_overall[9,3:]/num_avg_overall[9,3:], 'g', label = '1K wei-1')
plt.plot(x_array[3:], num_std_overall[10,3:]/num_avg_overall[10,3:], 'r', label = '1K wei-2')
plt.plot(x_array[3:], num_std_overall[11,3:]/num_avg_overall[11,3:], 'k', label = '1K wei-3')
plt.plot(x_array[3:], num_std_overall[12,3:]/num_avg_overall[12,3:], 'm', label = '1K wei-4')
plt.plot(x_array[3:], num_std_overall[13,3:]/num_avg_overall[13,3:], 'y', label = '1K mfa')
plt.xlabel("ensemble size ")
plt.ylabel("average std number concentration")
plt.grid()
plt.legend(loc = 'lower right')
fig = plt.gcf()
fig.savefig("figs/number_std_ensem_1K.pdf")

plt.clf()
plt.plot(x_array[1:], mass_avg_overall[7,1:], 'c', label = '1K wei+1')
plt.plot(x_array[1:], mass_avg_overall[8,1:], 'b', label = '1K flat')
plt.plot(x_array[1:], mass_avg_overall[9,1:], 'g', label = '1K wei-1')
plt.plot(x_array[1:], mass_avg_overall[10,1:], 'r', label = '1K wei-2')
plt.plot(x_array[1:], mass_avg_overall[11,1:], 'k', label = '1K wei-3')
plt.plot(x_array[1:], mass_avg_overall[12,1:], 'm', label = '1K wei-4')
plt.plot(x_array[1:], mass_avg_overall[13,1:], 'y', label = '1K mfa')
plt.xlabel("ensemble size ")
plt.ylabel("average mass concentration in \mu g m^{-3}")
plt.grid()
plt.legend(loc = 'lower right')
fig = plt.gcf()
fig.savefig("figs/mass_ensem_1K.pdf")

plt.clf()
plt.plot(x_array[3:], mass_std_overall[7,3:]/mass_avg_overall[7,3:], 'c', label = '1K wei+1')
plt.plot(x_array[3:], mass_std_overall[8,3:]/mass_avg_overall[8,3:], 'b', label = '1K flat')
plt.plot(x_array[3:], mass_std_overall[9,3:]/mass_avg_overall[9,3:], 'g', label = '1K wei-1')
plt.plot(x_array[3:], mass_std_overall[10,3:]/mass_avg_overall[10,3:], 'r', label = '1K wei-2')
plt.plot(x_array[3:], mass_std_overall[11,3:]/mass_avg_overall[11,3:], 'k', label = '1K wei-3')
plt.plot(x_array[3:], mass_std_overall[12,3:]/mass_avg_overall[12,3:], 'm', label = '1K wei-4')
plt.plot(x_array[3:], mass_std_overall[13,3:]/mass_avg_overall[13,3:], 'y', label = '1K mfa')
plt.grid()
plt.xlabel("ensemble size ")
plt.ylabel("average std mass concentration")
plt.legend(loc = 'lower right')
fig = plt.gcf()
fig.savefig("figs/mass_std_ensem_1K.pdf")


    



