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

#num_avg_100K_4 = np.loadtxt("data/average_number_100K_wei-4.txt")
#num_std_100K_4 = np.loadtxt("data/std_number_100K_wei-4.txt")
#mass_avg_100K_4 = np.loadtxt("data/average_mass_100K_wei-4.txt")
#mass_std_100K_4 = np.loadtxt("data/std_mass_100K_wei-4.txt")

num_avg_100K_3 = np.loadtxt("data/average_number_100K_wei-3.txt")
num_std_100K_3 = np.loadtxt("data/std_number_100K_wei-3.txt")
mass_avg_100K_3 = np.loadtxt("data/average_mass_100K_wei-3.txt")
mass_std_100K_3 = np.loadtxt("data/std_mass_100K_wei-3.txt")

num_avg_100K_2 = np.loadtxt("data/average_number_100K_wei-2.txt")
num_std_100K_2 = np.loadtxt("data/std_number_100K_wei-2.txt")
mass_avg_100K_2 = np.loadtxt("data/average_mass_100K_wei-2.txt")
mass_std_100K_2 = np.loadtxt("data/std_mass_100K_wei-2.txt")

num_avg_100K_1 = np.loadtxt("data/average_number_100K_wei-1.txt")
num_std_100K_1 = np.loadtxt("data/std_number_100K_wei-1.txt")
mass_avg_100K_1 = np.loadtxt("data/average_mass_100K_wei-1.txt")
mass_std_100K_1 = np.loadtxt("data/std_mass_100K_wei-1.txt")

num_avg_100K_flat = np.loadtxt("data/average_number_100K_flat.txt")
num_std_100K_flat = np.loadtxt("data/std_number_100K_flat.txt")
mass_avg_100K_flat = np.loadtxt("data/average_mass_100K_flat.txt")
mass_std_100K_flat = np.loadtxt("data/std_mass_100K_flat.txt")

num_avg_10K_4 = np.loadtxt("data/average_number_10K_wei-4.txt")
num_std_10K_4 = np.loadtxt("data/std_number_10K_wei-4.txt")
mass_avg_10K_4 = np.loadtxt("data/average_mass_10K_wei-4.txt")
mass_std_10K_4 = np.loadtxt("data/std_mass_10K_wei-4.txt")

num_avg_10K_3 = np.loadtxt("data/average_number_10K_wei-3.txt")
num_std_10K_3 = np.loadtxt("data/std_number_10K_wei-3.txt")
mass_avg_10K_3 = np.loadtxt("data/average_mass_10K_wei-3.txt")
mass_std_10K_3 = np.loadtxt("data/std_mass_10K_wei-3.txt")

num_avg_10K_2 = np.loadtxt("data/average_number_10K_wei-2.txt")
num_std_10K_2 = np.loadtxt("data/std_number_10K_wei-2.txt")
mass_avg_10K_2 = np.loadtxt("data/average_mass_10K_wei-2.txt")
mass_std_10K_2 = np.loadtxt("data/std_mass_10K_wei-2.txt")

num_avg_10K_1 = np.loadtxt("data/average_number_10K_wei-1.txt")
num_std_10K_1 = np.loadtxt("data/std_number_10K_wei-1.txt")
mass_avg_10K_1 = np.loadtxt("data/average_mass_10K_wei-1.txt")
mass_std_10K_1 = np.loadtxt("data/std_mass_10K_wei-1.txt")

num_avg_10K_flat = np.loadtxt("data/average_number_10K_flat.txt")
num_std_10K_flat = np.loadtxt("data/std_number_10K_flat.txt")
mass_avg_10K_flat = np.loadtxt("data/average_mass_10K_flat.txt")
mass_std_10K_flat = np.loadtxt("data/std_mass_10K_flat.txt")

# average std over time
array_num_std_avg_100K = np.zeros([5])
array_num_std_avg_10K = np.zeros([5])
array_mass_std_avg_100K = np.zeros([5])
array_mass_std_avg_10K = np.zeros([5])

array_num_std_avg_100K[0] = np.average(num_std_100K_flat/num_avg_100K_flat)
array_num_std_avg_100K[1] = np.average(num_std_100K_1/num_avg_100K_1)
array_num_std_avg_100K[2] = np.average(num_std_100K_2/num_avg_100K_2)
array_num_std_avg_100K[3] = np.average(num_std_100K_3/num_avg_100K_3)
#array_num_std_avg_100K[4] = np.average(num_std_100K_4/num_avg_100K_4)

array_num_std_avg_10K[0] = np.average(num_std_10K_flat/num_avg_10K_flat)
array_num_std_avg_10K[1] = np.average(num_std_10K_1/num_avg_10K_1)
array_num_std_avg_10K[2] = np.average(num_std_10K_2/num_avg_10K_2)
array_num_std_avg_10K[3] = np.average(num_std_10K_3/num_avg_10K_3)
array_num_std_avg_10K[4] = np.average(num_std_10K_4/num_avg_10K_4)

array_mass_std_avg_100K[0] = np.average(mass_std_100K_flat/mass_avg_100K_flat)
array_mass_std_avg_100K[1] = np.average(mass_std_100K_1/mass_avg_100K_1)
array_mass_std_avg_100K[2] = np.average(mass_std_100K_2/mass_avg_100K_2)
array_mass_std_avg_100K[3] = np.average(mass_std_100K_3/mass_avg_100K_3)
#array_mass_std_avg_100K[4] = np.average(mass_std_100K_4/mass_avg_100K_4)

array_mass_std_avg_10K[0] = np.average(mass_std_10K_flat/mass_avg_10K_flat)
array_mass_std_avg_10K[1] = np.average(mass_std_10K_1/mass_avg_10K_1)
array_mass_std_avg_10K[2] = np.average(mass_std_10K_2/mass_avg_10K_2)
array_mass_std_avg_10K[3] = np.average(mass_std_10K_3/mass_avg_10K_3)
array_mass_std_avg_10K[4] = np.average(mass_std_10K_4/mass_avg_10K_4)

x_array = [0, -1, -2, -3, -4]
print x_array
print array_num_std_avg_100K, array_num_std_avg_10K, array_mass_std_avg_100K, array_mass_std_avg_10K
plt.clf()
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
plt.plot(x_array, array_num_std_avg_100K, 'r-', label = 'number 100K')
plt.plot(x_array, array_num_std_avg_10K, 'b-', label = 'number 10K')
plt.plot(x_array, array_mass_std_avg_100K, 'r--', label = 'mass 100K')
plt.plot(x_array, array_mass_std_avg_10K, 'b--', label = 'mass 10K')
plt.legend(loc = 'upper right')
plt.xlabel("xxx ")
plt.ylabel("yyy")
plt.title("zzz")
fig = plt.gcf()
fig.savefig("figs/std_total_averages2.pdf")



plt.clf()
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
plt.plot(time_array, mass_std_100K_flat/mass_avg_100K_flat, '--', label = '100K 0')
plt.plot(time_array, mass_std_100K_1/mass_avg_100K_1, '--', label = '100K -1')
plt.plot(time_array, mass_std_100K_2/mass_avg_100K_2, '--', label = '100K -2')
plt.plot(time_array, mass_std_100K_3/mass_avg_100K_3, '--', label = '100K -3')
plt.plot(time_array, mass_std_10K_flat/mass_avg_10K_flat, label = '10K 0')
plt.plot(time_array, mass_std_10K_1/mass_avg_10K_1, label = '10K -1')
plt.plot(time_array, mass_std_10K_2/mass_avg_10K_2, label = '10K -2')
plt.plot(time_array, mass_std_10K_3/mass_avg_10K_3, label = '10K -3')
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
plt.plot(time_array, num_std_100K_flat/num_avg_100K_flat, '--', label = '100K 0')
plt.plot(time_array, num_std_100K_1/num_avg_100K_1, '--', label = '100K -1')
plt.plot(time_array, num_std_100K_2/num_avg_100K_2, '--', label = '100K -2')
plt.plot(time_array, num_std_100K_3/num_avg_100K_3, '--', label = '100K -3')
plt.plot(time_array, num_std_10K_flat/num_avg_10K_flat, label = '10K 0')
plt.plot(time_array, num_std_10K_1/num_avg_10K_1, label = '10K -1')
plt.plot(time_array, num_std_10K_2/num_avg_10K_2, label = '10K -2')
plt.plot(time_array, num_std_10K_3/num_avg_10K_3, label = '10K -3')
plt.legend(loc = 'upper right')
plt.xlabel("time ")
plt.ylabel("std/avg")
plt.title("Number")
fig = plt.gcf()
fig.savefig("figs/std_total_number.pdf")




    



