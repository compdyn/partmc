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

time_array = np.loadtxt("data/time_array_100K_wei-3.txt")
i_counter = 0
i_counter_p = 0

array_num_std_avg = np.zeros([3,5])
array_mass_std_avg = np.zeros([3,5])

num_avg = np.zeros([3,5,25])
num_std = np.zeros([3,5,25])
mass_avg = np.zeros([3,5,25])
mass_std = np.zeros([3,5,25])

for counter_p in ["1K", "10K", "100K"]:
    print counter_p, i_counter_p
    i_counter = 0
    for counter in ["flat", "wei-1", "wei-2", "wei-3", "wei-4"]:
        print counter, i_counter
        f1 = "data/average_number_%s_%s.txt" % (counter_p, counter)
        f2 = "data/std_number_%s_%s.txt" % (counter_p, counter)
        f3 = "data/average_mass_%s_%s.txt" % (counter_p, counter)
        f4 = "data/std_mass_%s_%s.txt" % (counter_p, counter)

        print f1, f2, f3, f4
        num_avg[i_counter_p,i_counter,:] = np.loadtxt(f1) 
        num_std[i_counter_p,i_counter,:] = np.loadtxt(f2) 
        mass_avg[i_counter_p,i_counter,:] = np.loadtxt(f3) 
        mass_std[i_counter_p,i_counter,:] = np.loadtxt(f4) 

# average std over time

        array_num_std_avg[i_counter_p,i_counter] = np.average(num_std[i_counter_p,i_counter,:]/num_avg[i_counter_p,i_counter,:])
        array_mass_std_avg[i_counter_p,i_counter] = np.average(mass_std[i_counter_p, i_counter,:]/mass_avg[i_counter_p,i_counter,:])
        i_counter = i_counter + 1
    i_counter_p = i_counter_p + 1

print 'array_num_std_avg ', array_num_std_avg

x_array = [0, -1, -2, -3, -4]

plt.clf()
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
plt.plot(x_array, array_num_std_avg[0,:], 'r-', label = 'number 1K')
plt.plot(x_array, array_num_std_avg[1,:], 'g-', label = 'number 10K')
plt.plot(x_array, array_num_std_avg[2,:], 'b-', label = 'number 100K')
plt.plot(x_array, array_mass_std_avg[0,:], 'r--', label = 'mass 1K')
plt.plot(x_array, array_mass_std_avg[1,:], 'g--', label = 'mass 10K')
plt.plot(x_array, array_mass_std_avg[2,:], 'b--', label = 'mass 100K')
plt.legend(loc = 'upper right')
plt.xlabel("weighting exponent ")
plt.ylabel("averaged std/avg")
fig = plt.gcf()
fig.savefig("figs/std_total_averages.pdf")

plt.clf()
a = plt.gca()
a.set_xscale("log")
a.set_yscale("log")
plt.plot(array_num_std_avg[:,0], array_mass_std_avg[:,0], 'xr-', label = '0')
plt.plot(array_num_std_avg[:,1], array_mass_std_avg[:,1], 'xg-', label = '-1')
plt.plot(array_num_std_avg[:,2], array_mass_std_avg[:,2], 'xb-', label = '-2')
plt.plot(array_num_std_avg[:,3], array_mass_std_avg[:,3], 'xm-', label = '-3')
plt.plot(array_num_std_avg[:,4], array_mass_std_avg[:,4], 'xk-', label = '-4')
plt.legend(loc = 'upper right')
plt.xlabel("number std/avg ")
plt.ylabel("mass std/avg")
fig = plt.gcf()
fig.savefig("figs/std_total_averages_mvn_1.pdf")

plt.clf()
a = plt.gca()
a.set_xscale("log")
a.set_yscale("log")
plt.plot(array_num_std_avg[0,:], array_mass_std_avg[0,:], 'xr-', label = '1K')
plt.plot(array_num_std_avg[1,:], array_mass_std_avg[1,:], 'xg-', label = '10K')
plt.plot(array_num_std_avg[2,:], array_mass_std_avg[2,:], 'xb-', label = '100K')

plt.legend(loc = 'upper right')
plt.xlabel("number std/avg ")
plt.ylabel("mass std/avg")
fig = plt.gcf()
fig.savefig("figs/std_total_averages_mvn_2.pdf")

plt.clf()
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
plt.plot(time_array, mass_std[0,0,:]/mass_avg[0,0,:], 'r-', label = '1K 0')
plt.plot(time_array, mass_std[0,1,:]/mass_avg[0,1,:], 'r-', label = '1K -1')
plt.plot(time_array, mass_std[0,2,:]/mass_avg[0,2,:], 'r-', label = '1K -2')
plt.plot(time_array, mass_std[0,3,:]/mass_avg[0,3,:], 'r-', label = '1K -3')
plt.plot(time_array, mass_std[1,0,:]/mass_avg[1,0,:], 'g-', label = '10K 0')
plt.plot(time_array, mass_std[1,1,:]/mass_avg[1,1,:], 'g-', label = '10K -1')
plt.plot(time_array, mass_std[1,2,:]/mass_avg[1,2,:], 'g-', label = '10K -2')
plt.plot(time_array, mass_std[1,3,:]/mass_avg[1,3,:], 'g-', label = '10K -3')
plt.plot(time_array, mass_std[2,0,:]/mass_avg[2,0,:], 'b-', label = '100K 0')
plt.plot(time_array, mass_std[2,1,:]/mass_avg[2,1,:], 'b-', label = '100K -1')
plt.plot(time_array, mass_std[2,2,:]/mass_avg[2,2,:], 'b-', label = '100K -2')
plt.plot(time_array, mass_std[2,3,:]/mass_avg[2,3,:], 'b-', label = '100K -3')
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
plt.plot(time_array, num_std[0,0,:]/num_avg[0,0,:], 'r-', label = '1K 0')
plt.plot(time_array, num_std[0,1,:]/num_avg[0,1,:], 'r-', label = '1K -1')
plt.plot(time_array, num_std[0,2,:]/num_avg[0,2,:], 'r-', label = '1K -2')
plt.plot(time_array, num_std[0,3,:]/num_avg[0,3,:], 'r-', label = '1K -3')
plt.plot(time_array, num_std[1,0,:]/num_avg[1,0,:], 'g-', label = '10K 0')
plt.plot(time_array, num_std[1,1,:]/num_avg[1,1,:], 'g-', label = '10K -1')
plt.plot(time_array, num_std[1,2,:]/num_avg[1,2,:], 'g-', label = '10K -2')
plt.plot(time_array, num_std[1,3,:]/num_avg[1,3,:], 'g-', label = '10K -3')
plt.plot(time_array, num_std[2,0,:]/num_avg[2,0,:], 'b-', label = '100K 0')
plt.plot(time_array, num_std[2,1,:]/num_avg[2,1,:], 'b-', label = '100K -1')
plt.plot(time_array, num_std[2,2,:]/num_avg[2,2,:], 'b-', label = '100K -2')
plt.plot(time_array, num_std[2,3,:]/num_avg[2,3,:], 'b-', label = '100K -3')
plt.legend(loc = 'upper right')
plt.xlabel("time ")
plt.ylabel("std/avg")
plt.title("Number")
fig = plt.gcf()
fig.savefig("figs/std_total_number.pdf")




    



