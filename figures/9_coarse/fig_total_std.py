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

array_num_std_avg = np.zeros([3,7]) # 3 for particle number cases, 7 for weighting cases
array_mass_std_avg = np.zeros([3,7])

num_avg = np.zeros([3,7,25])
num_std = np.zeros([3,7,25])
mass_avg = np.zeros([3,7,25])
mass_std = np.zeros([3,7,25])

for counter_p in ["1K", "10K", "100K"]:
    print counter_p, i_counter_p
    i_counter = 0
    for counter in ["wei\+1", "flat", "wei-1", "wei-2", "wei-3", "wei-4", "mfa"]:
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

x_array = [1, 0, -1, -2, -3, -4]

print 'test ',array_num_std_avg[0,0:5], array_num_std_avg[0,6]
plt.clf()
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
plt.plot(x_array, array_num_std_avg[0,0:6], 'r-x', label = 'number 1K')
plt.plot(-3, array_num_std_avg[0,6], 'ro', label = 'number 1K mfa')
plt.plot(x_array, array_num_std_avg[1,0:6], 'g-x', label = 'number 10K')
plt.plot(-3, array_num_std_avg[1,6], 'go', label = 'number 10K mfa')
plt.plot(x_array, array_num_std_avg[2,0:6], 'b-x', label = 'number 100K')
plt.plot(-3, array_num_std_avg[2,6], 'bo', label = 'number 100K mfa')
plt.plot(x_array, array_mass_std_avg[0,0:6], 'r--x', label = 'mass 1K')
plt.plot(-3, array_mass_std_avg[0,6], 'ro', label = 'mass 1K mfa')
plt.plot(x_array, array_mass_std_avg[1,0:6], 'g--x', label = 'mass 10K')
plt.plot(-3, array_mass_std_avg[1,6], 'go', label = 'mass 10K mfa')
plt.plot(x_array, array_mass_std_avg[2,0:6], 'b--x', label = 'mass 100K')
plt.plot(-3, array_mass_std_avg[2,6], 'bo', label = 'mass 100K mfa')
plt.grid()
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
#plt.plot(array_num_std_avg[:,2], array_mass_std_avg[:,2], 'xb-', label = '-2')
plt.plot(array_num_std_avg[:,3], array_mass_std_avg[:,3], 'xm-', label = '-3')
#plt.plot(array_num_std_avg[:,4], array_mass_std_avg[:,4], 'xk-', label = '-4')
plt.grid()
plt.legend(loc = 'upper right')
plt.xlabel("number std/avg ")
plt.ylabel("mass std/avg")
fig = plt.gcf()
fig.savefig("figs/std_total_averages_mvn_1.pdf")

plt.clf()
a = plt.gca()
a.set_xscale("log")
a.set_yscale("log")
plt.plot(array_num_std_avg[0,0:6], array_mass_std_avg[0,0:6], 'xr-', label = '1K')
plt.plot(array_num_std_avg[0,6], array_mass_std_avg[0,6], 'or', label = '1K mfa')
plt.plot(array_num_std_avg[1,0:6], array_mass_std_avg[1,0:6], 'xg-', label = '10K')
plt.plot(array_num_std_avg[1,6], array_mass_std_avg[1,6], 'og', label = '10K mfa')
plt.plot(array_num_std_avg[2,0:6], array_mass_std_avg[2,0:6], 'xb-', label = '100K')
plt.plot(array_num_std_avg[2,6], array_mass_std_avg[2,6], 'ob', label = '100K mfa')
plt.grid()
plt.legend(loc = 'upper right')
plt.xlabel("number std/avg ")
plt.ylabel("mass std/avg")
fig = plt.gcf()
fig.savefig("figs/std_total_averages_mvn_2.pdf")

plt.clf()
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
plt.plot(time_array, mass_std[0,0,:]/mass_avg[0,0,:], 'r-', label = '1K +1')
plt.plot(time_array, mass_std[0,1,:]/mass_avg[0,1,:], 'b-', label = '1K 0')
plt.plot(time_array, mass_std[0,2,:]/mass_avg[0,2,:], 'g-', label = '1K -1')
plt.plot(time_array, mass_std[0,3,:]/mass_avg[0,3,:], 'm-', label = '1K -2')
plt.plot(time_array, mass_std[0,4,:]/mass_avg[0,4,:], 'y-', label = '1K -3')
plt.plot(time_array, mass_std[0,5,:]/mass_avg[0,5,:], 'k-', label = '1K -4')
plt.plot(time_array, mass_std[0,6,:]/mass_avg[0,6,:], 'y-.', label = '1K mfa', linewidth=4)
plt.plot(time_array, mass_std[1,0,:]/mass_avg[1,0,:], 'r--', label = '10K +1')
plt.plot(time_array, mass_std[1,1,:]/mass_avg[1,1,:], 'b--', label = '10K 0')
plt.plot(time_array, mass_std[1,2,:]/mass_avg[1,2,:], 'g--', label = '10K -1')
plt.plot(time_array, mass_std[1,3,:]/mass_avg[1,3,:], 'm--', label = '10K -2')
plt.plot(time_array, mass_std[1,4,:]/mass_avg[1,4,:], 'y--', label = '10K -3')
plt.plot(time_array, mass_std[1,5,:]/mass_avg[1,5,:], 'k--', label = '10K -4')
plt.plot(time_array, mass_std[1,6,:]/mass_avg[1,6,:], 'y:', label = '10K mfa', linewidth=4)
plt.plot(time_array, mass_std[2,0,:]/mass_avg[2,0,:], 'r-', label = '100K +1')
plt.plot(time_array, mass_std[2,1,:]/mass_avg[2,1,:], 'b-', label = '100K 0')
plt.plot(time_array, mass_std[2,2,:]/mass_avg[2,2,:], 'g-', label = '100K -1')
plt.plot(time_array, mass_std[2,3,:]/mass_avg[2,3,:], 'm-', label = '100K -2')
plt.plot(time_array, mass_std[2,4,:]/mass_avg[2,4,:], 'y-', label = '100K -3')
plt.plot(time_array, mass_std[2,5,:]/mass_avg[2,5,:], 'k-', label = '100K -4')
plt.plot(time_array, mass_std[2,6,:]/mass_avg[2,6,:], 'y-.', label = '100K mfa', linewidth = 4)
plt.grid()
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
plt.plot(time_array, num_std[0,0,:]/num_avg[0,0,:], 'r-', label = '1K +1')
plt.plot(time_array, num_std[0,1,:]/num_avg[0,1,:], 'b-', label = '1K 0')
plt.plot(time_array, num_std[0,2,:]/num_avg[0,2,:], 'g-', label = '1K -1')
plt.plot(time_array, num_std[0,3,:]/num_avg[0,3,:], 'm-', label = '1K -2')
plt.plot(time_array, num_std[0,4,:]/num_avg[0,4,:], 'y-', label = '1K -3')
plt.plot(time_array, num_std[0,5,:]/num_avg[0,5,:], 'k-', label = '1K -4')
plt.plot(time_array, num_std[0,6,:]/num_avg[0,6,:], 'y-.', label = '1K mfa', linewidth=4)
plt.plot(time_array, num_std[1,0,:]/num_avg[1,0,:], 'r--', label = '10K +1')
plt.plot(time_array, num_std[1,1,:]/num_avg[1,1,:], 'b--', label = '10K 0')
plt.plot(time_array, num_std[1,2,:]/num_avg[1,2,:], 'g--', label = '10K -1')
plt.plot(time_array, num_std[1,3,:]/num_avg[1,3,:], 'm--', label = '10K -2')
plt.plot(time_array, num_std[1,4,:]/num_avg[1,4,:], 'y--', label = '10K -3')
plt.plot(time_array, num_std[1,5,:]/num_avg[1,5,:], 'k--', label = '10K -4')
plt.plot(time_array, num_std[1,6,:]/num_avg[1,6,:], 'y:', label = '10K mfa', linewidth=4)
plt.plot(time_array, num_std[2,0,:]/num_avg[2,0,:], 'r-', label = '100K +1')
plt.plot(time_array, num_std[2,1,:]/num_avg[2,1,:], 'b-', label = '100K 0')
plt.plot(time_array, num_std[2,2,:]/num_avg[2,2,:], 'g-', label = '100K -1')
plt.plot(time_array, num_std[2,3,:]/num_avg[2,3,:], 'm-', label = '100K -2')
plt.plot(time_array, num_std[2,4,:]/num_avg[2,4,:], 'y-', label = '100K -3')
plt.plot(time_array, num_std[2,5,:]/num_avg[2,5,:], 'k-', label = '100K -4')
plt.plot(time_array, num_std[2,6,:]/num_avg[2,6,:], 'y-.', label = '100K mfa', linewidth = 4)
plt.grid()
plt.legend(loc = 'upper right')
plt.xlabel("time ")
plt.ylabel("std/avg")
plt.title("Number")
fig = plt.gcf()
fig.savefig("figs/std_total_number.pdf")




    



