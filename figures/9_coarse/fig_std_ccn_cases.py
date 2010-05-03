#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
import config

ccn_std_ss1 = np.loadtxt("data/ccn_std_ss1.txt")
ccn_std_ss2 = np.loadtxt("data/ccn_std_ss2.txt")
ccn_std_ss3 = np.loadtxt("data/ccn_std_ss3.txt")
ccn_std_ss4 = np.loadtxt("data/ccn_std_ss4.txt")

ccn_average_ss1 = np.loadtxt("data/ccn_average_ss1.txt")
ccn_average_ss2 = np.loadtxt("data/ccn_average_ss2.txt")
ccn_average_ss3 = np.loadtxt("data/ccn_average_ss3.txt")
ccn_average_ss4 = np.loadtxt("data/ccn_average_ss4.txt")

ccn_array_ss1 = np.zeros([25,5])
ccn_array_ss2 = np.zeros([25,5])
ccn_array_ss3 = np.zeros([25,5])
ccn_array_ss4 = np.zeros([25,5])
for i in range(0,5):
    ccn_array_ss1[:,i] = ccn_std_ss1[:,i] / ccn_average_ss1[:,i]
    ccn_array_ss2[:,i] = ccn_std_ss2[:,i] / ccn_average_ss2[:,i]
    ccn_array_ss3[:,i] = ccn_std_ss3[:,i] / ccn_average_ss3[:,i]
    ccn_array_ss4[:,i] = ccn_std_ss4[:,i] / ccn_average_ss4[:,i]

print ccn_array_ss1
plt.clf()
plt.plot(ccn_array_ss1[:,0], 'b-', label = '10K flat')
plt.plot(ccn_array_ss1[:,1], 'g-', label = '10K -1')
plt.plot(ccn_array_ss1[:,2], 'r-', label = '10K -2')
plt.plot(ccn_array_ss1[:,3], 'b--', label = '10K -3')
plt.plot(ccn_array_ss1[:,4], 'g--', label = '10K -4')
#plt.axis([0, 24, 0, 1])
plt.legend(loc = 'upper right')
plt.xlabel("time (hours)")
plt.ylabel("std/avg")
plt.title('S_c = 0.001')
fig = plt.gcf()
fig.savefig("figs/ccn_std_10K_ss1.pdf")

plt.clf()
plt.plot(ccn_array_ss2[:,0], 'b-', label = '10K flat')
plt.plot(ccn_array_ss2[:,1], 'g-', label = '10K -1')
plt.plot(ccn_array_ss2[:,2], 'r-', label = '10K -2')
plt.plot(ccn_array_ss2[:,3], 'b--', label = '10K -3')
plt.plot(ccn_array_ss2[:,4], 'g--', label = '10K -4')
#plt.axis([0, 24, 0, 1])
plt.legend(loc = 'upper right')
plt.xlabel("time (hours)")
plt.ylabel("std/avg")
plt.title('S_c = 0.01')
fig = plt.gcf()
fig.savefig("figs/ccn_std_10K_ss2.pdf")

plt.clf()
plt.plot(ccn_array_ss3[:,0], 'b-', label = '10K flat')
plt.plot(ccn_array_ss3[:,1], 'g-', label = '10K -1')
plt.plot(ccn_array_ss3[:,2], 'r-', label = '10K -2')
plt.plot(ccn_array_ss3[:,3], 'b--', label = '10K -3')
plt.plot(ccn_array_ss3[:,4], 'g--', label = '10K -4')
#plt.axis([0, 24, 0, 1])
plt.legend(loc = 'upper right')
plt.xlabel("time (hours)")
plt.ylabel("std/avg")
plt.title('S_c = 0.1')
fig = plt.gcf()
fig.savefig("figs/ccn_std_10K_ss3.pdf")

plt.clf()
plt.plot(ccn_array_ss4[:,0], 'b-', label = '10K flat')
plt.plot(ccn_array_ss4[:,1], 'g-', label = '10K -1')
plt.plot(ccn_array_ss4[:,2], 'r-', label = '10K -2')
plt.plot(ccn_array_ss4[:,3], 'b--', label = '10K -3')
plt.plot(ccn_array_ss4[:,4], 'g--', label = '10K -4')
#plt.axis([0, 24, 0, 1])
plt.legend(loc = 'upper right')
plt.xlabel("time (hours)")
plt.ylabel("std/avg")
plt.title('S_c = 0.5')
fig = plt.gcf()
fig.savefig("figs/ccn_std_10K_ss4.pdf")

