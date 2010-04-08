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
i_counter = 0

array_num_std_avg_100K = np.zeros([5])
array_num_std_avg_10K = np.zeros([5])
array_mass_std_avg_100K = np.zeros([5])
array_mass_std_avg_10K = np.zeros([5])

num_avg_10K = np.zeros([5,25])
num_std_10K = np.zeros([5,25])
mass_avg_10K = np.zeros([5,25])
mass_std_10K = np.zeros([5,25])

num_avg_100K = np.zeros([5,25])
num_std_100K = np.zeros([5,25])
mass_avg_100K = np.zeros([5,25])
mass_std_100K = np.zeros([5,25])

for counter in ["flat", "wei-1", "wei-2", "wei-3", "wei-4"]:
    f1 = "data/average_number_10K_%s.txt" % counter
    f2 = "data/std_number_10K_%s.txt" % counter
    f3 = "data/average_mass_10K_%s.txt" % counter
    f4 = "data/std_mass_10K_%s.txt" % counter

    num_avg_10K[i_counter,:] = np.loadtxt(f1) 
    num_std_10K[i_counter,:] = np.loadtxt(f2) 
    mass_avg_10K[i_counter,:] = np.loadtxt(f3) 
    mass_std_10K[i_counter,:] = np.loadtxt(f4) 

    f5 = "data/average_number_100K_%s.txt" % counter
    f6 = "data/std_number_100K_%s.txt" % counter
    f7 = "data/average_mass_100K_%s.txt" % counter
    f8 = "data/std_mass_100K_%s.txt" % counter

    num_avg_100K[i_counter,:] = np.loadtxt(f5) 
    num_std_100K[i_counter,:] = np.loadtxt(f6) 
    mass_avg_100K[i_counter,:] = np.loadtxt(f7) 
    mass_std_100K[i_counter,:] = np.loadtxt(f8) 

# average std over time

    array_num_std_avg_10K[i_counter] = np.average(num_std_10K[i_counter,:]/num_avg_10K[i_counter,:])
    array_mass_std_avg_10K[i_counter] = np.average(mass_std_10K[i_counter,:]/mass_avg_10K[i_counter,:])
    array_num_std_avg_100K[i_counter] = np.average(num_std_100K[i_counter,:]/num_avg_100K[i_counter,:])
    array_mass_std_avg_100K[i_counter] = np.average(mass_std_100K[i_counter,:]/mass_avg_100K[i_counter,:])
    i_counter = i_counter + 1

x_array = [0, -1, -2, -3, -4]
print 'array_num_std_avg_10K ', array_num_std_avg_10K

plt.clf()
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
plt.plot(x_array, array_num_std_avg_100K, 'r-', label = 'number 100K')
plt.plot(x_array, array_num_std_avg_10K, 'b-', label = 'number 10K')
plt.plot(x_array, array_mass_std_avg_100K, 'r--', label = 'mass 100K')
plt.plot(x_array, array_mass_std_avg_10K, 'b--', label = 'mass 10K')
plt.legend(loc = 'upper right')
plt.xlabel("weighting exponent ")
plt.ylabel("averaged std/avg")
fig = plt.gcf()
fig.savefig("figs/std_total_averages.pdf")

plt.clf()
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
plt.plot(time_array, mass_std_100K[0,:]/mass_avg_100K[0,:], '--', label = '100K 0')
plt.plot(time_array, mass_std_100K[1,:]/mass_avg_100K[1,:], '--', label = '100K -1')
plt.plot(time_array, mass_std_100K[2,:]/mass_avg_100K[2,:], '--', label = '100K -2')
plt.plot(time_array, mass_std_100K[3,:]/mass_avg_100K[3,:], '--', label = '100K -3')
plt.plot(time_array, mass_std_10K[0,:]/mass_avg_10K[0,:], label = '10K 0')
plt.plot(time_array, mass_std_10K[1,:]/mass_avg_10K[1,:], label = '10K -1')
plt.plot(time_array, mass_std_10K[2,:]/mass_avg_10K[2,:], label = '10K -2')
plt.plot(time_array, mass_std_10K[3,:]/mass_avg_10K[3,:], label = '10K -3')
plt.legend(loc = 'upper right')
plt.xlabel("time ")
plt.ylabel("std/avg")
plt.title("Mass")
fig = plt.gcf()
fig.savefig("figs/std_total_mass.pdf")

plt.clf()
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
plt.plot(time_array, num_std_100K[0,:]/num_avg_100K[0,:], '--', label = '100K 0')
plt.plot(time_array, num_std_100K[1,:]/num_avg_100K[1,:], '--', label = '100K -1')
plt.plot(time_array, num_std_100K[2,:]/num_avg_100K[2,:], '--', label = '100K -2')
plt.plot(time_array, num_std_100K[3,:]/num_avg_100K[3,:], '--', label = '100K -3')
plt.plot(time_array, num_std_10K[0,:]/num_avg_10K[0,:], label = '10K 0')
plt.plot(time_array, num_std_10K[1,:]/num_avg_10K[1,:], label = '10K -1')
plt.plot(time_array, num_std_10K[2,:]/num_avg_10K[2,:], label = '10K -2')
plt.plot(time_array, num_std_10K[3,:]/num_avg_10K[3,:], label = '10K -3')
plt.legend(loc = 'upper right')
plt.xlabel("time ")
plt.ylabel("std/avg")
plt.title("Number")
fig = plt.gcf()
fig.savefig("figs/std_total_number.pdf")




    



