#!/usr/bin/env python

import scipy.io
import scipy.stats
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc
import config

#calculation of confidence interval

c = config.c_value
n = config.i_loop_max
r = scipy.stats.t.ppf((1 + c) / 2, n - 1)
conf_factor = r / np.sqrt(n)

print 'c, n, r', c, n, r, r/np.sqrt(n)

num_avg_overall = np.zeros([21,config.i_loop_max])
mass_avg_overall = np.zeros([21,config.i_loop_max])
num_std_overall = np.zeros([21,config.i_loop_max])
mass_std_overall = np.zeros([21,config.i_loop_max])

error_ratio_1K_num = np.zeros([7])
error_ratio_1K_mass = np.zeros([7])

error_ratio_10K_num = np.zeros([7])
error_ratio_10K_mass = np.zeros([7])

i_counter = 0
for counter in ["1K_wei\\+1","10K_wei\\+1", "100K_wei\\+1",
		"1K_flat", "10K_flat", "100K_flat", 
		"1K_wei-1", "10K_wei-1", "100K_wei-1", 
		"1K_wei-2", "10K_wei-2", "100K_wei-2", 
		"1K_wei-3", "10K_wei-3", "100K_wei-3", 
		"1K_wei-4", "10K_wei-4", "100K_wei-4", 
		"1K_mfa", "10K_mfa", "100K_mfa", ]:
    f1 = "data/ensemble_number_%s.txt" % counter
    f2 = "data/ensemble_mass_%s.txt" % counter
    f3 = "data/ensemble_number_std_%s.txt" % counter
    f4 = "data/ensemble_mass_std_%s.txt" % counter

    num_avg_overall[i_counter,:] = np.loadtxt(f1)
    mass_avg_overall[i_counter,:] = np.loadtxt(f2)
    num_std_overall[i_counter,:] = np.loadtxt(f3)
    mass_std_overall[i_counter,:] = np.loadtxt(f4)

    i_counter +=1

x_array = [1000, 10000, 1e5]

print 'check ', x_array[0],x_array[1]

i_weight = 0
for weight in ["wei+1", "flat", "wei-1", "wei-2", "wei-3", "wei-4", "mfa"]:
    print "weight ", i_weight, weight
    print "index ", i_weight*3, i_weight*3+1, i_weight*3 + 2
    plt.clf()
    plt.errorbar(x_array[0], num_avg_overall[i_weight*3,99], r*num_std_overall[i_weight*3,99], marker='s', mfc='b', fmt='o', ecolor='g')
    plt.errorbar(x_array[0], num_avg_overall[i_weight*3,99], conf_factor * num_std_overall[i_weight*3,99], 
                 marker='s', mfc='b', fmt='o', ecolor='r')
    plt.errorbar(x_array[1], num_avg_overall[1+i_weight*3,99], r*num_std_overall[1+i_weight*3,99], marker='s', mfc='b',fmt='o', ecolor='g')
    plt.errorbar(x_array[1], num_avg_overall[1+i_weight*3,99], conf_factor * num_std_overall[1+i_weight*3,99], 
                 marker='s', mfc='b',fmt='o', ecolor='r')
    plt.errorbar(x_array[2], num_avg_overall[2+i_weight*3,99], r*num_std_overall[2+i_weight*3,99], marker='s', mfc='b',fmt='o', ecolor='g')
    plt.errorbar(x_array[2], num_avg_overall[2+i_weight*3,99], conf_factor * num_std_overall[2+i_weight*3,99], 
                 marker='s', mfc='b',fmt='o', ecolor='r')
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("linear")
    plt.xlim([500, 2e5])
    plt.grid(True)
    plt.xlabel("particle number ")
    plt.ylabel("mean number in m^{-3}")
    fig = plt.gcf()
    out_file = "figs/num_mean_v_particlen_"+weight+".pdf"
    fig.savefig(out_file)
    i_weight += 1

i_weight = 0
for weight in ["wei+1", "flat", "wei-1", "wei-2", "wei-3", "wei-4", "mfa"]:
    plt.clf()
    plt.errorbar(x_array[0], mass_avg_overall[i_weight*3,99], r*mass_std_overall[i_weight*3,99], marker='s', mfc='b', fmt='o', ecolor='g')
    plt.errorbar(x_array[0], mass_avg_overall[i_weight*3,99], conf_factor * mass_std_overall[i_weight*3,99],
                 marker='s', mfc='b', fmt='o', ecolor='r')
    plt.errorbar(x_array[1], mass_avg_overall[1+i_weight*3,99], r*mass_std_overall[1+i_weight*3,99], marker='s', mfc='b',fmt='o', ecolor='g')
    plt.errorbar(x_array[1], mass_avg_overall[1+i_weight*3,99], conf_factor * mass_std_overall[1+i_weight*3,99],
                 marker='s', mfc='b',fmt='o', ecolor='r')
    plt.errorbar(x_array[2], mass_avg_overall[2+i_weight*3,99], r*mass_std_overall[2+i_weight*3,99], marker='s', mfc='b',fmt='o', ecolor='g')
    plt.errorbar(x_array[2], mass_avg_overall[2+i_weight*3,99], conf_factor * mass_std_overall[2+i_weight*3,99],
                 marker='s', mfc='b',fmt='o', ecolor='r')
    a = plt.gca()
    a.set_xscale("log")
    a.set_yscale("linear")
    plt.xlim([500, 2e5])
    plt.grid(True)
    plt.xlabel("particle number ")
    plt.ylabel("mean mass in kg m^{-3}")
    fig = plt.gcf()
    out_file = "figs/mass_mean_v_particlen_"+weight+".pdf"
    fig.savefig(out_file)
    i_weight += 1



i_weight = 0
for weight in ["wei+1", "flat", "wei-1", "wei-2", "wei-3", "wei-4", "mfa"]:
    index_1K = 3*i_weight
    index_10K = 3*i_weight + 1
    print "ER ", weight, index_1K + 1, index_1K, index_10K+1, index_10K
    print "values ", mass_avg_overall[index_1K+1, 99], mass_avg_overall[index_1K, 99], mass_std_overall[index_1K,99]
    print "values ", mass_avg_overall[index_10K+1, 99], mass_avg_overall[index_10K, 99], mass_std_overall[index_10K,99]

    error_ratio_1K_num[i_weight] = (num_avg_overall[index_1K+1, 99] - num_avg_overall[index_1K, 99])/(r*num_std_overall[index_1K,99])    
    error_ratio_10K_num[i_weight] = (num_avg_overall[index_10K+1, 99] - num_avg_overall[index_10K, 99])/(r*num_std_overall[index_10K,99]) 

    error_ratio_1K_mass[i_weight] = (mass_avg_overall[index_1K+1, 99] - mass_avg_overall[index_1K, 99])/(r*mass_std_overall[index_1K,99])       
    error_ratio_10K_mass[i_weight] = (mass_avg_overall[index_10K+1, 99] - mass_avg_overall[index_10K, 99])/(r*mass_std_overall[index_10K,99])
    i_weight += 1

weight_array = [1, 0, -1, -2, -3, -4]
print"error_ratio ", weight_array, error_ratio_1K_num[0:7], error_ratio_10K_num
plt.clf()
plt.plot(weight_array, error_ratio_1K_num[0:6], label = "num, 1K")
plt.plot(weight_array, error_ratio_10K_num[0:6], label = "num, 10K")
plt.plot(weight_array, error_ratio_1K_mass[0:6], label = "mass, 1K")
plt.plot(weight_array, error_ratio_10K_mass[0:6], label = "mass, 10K")
plt.legend()
plt.grid(True)
plt.xlabel("\alpha ")
plt.ylabel("ER")
fig = plt.gcf()
fig.savefig("error_ratio.pdf")








