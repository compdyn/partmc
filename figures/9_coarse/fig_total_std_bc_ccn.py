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

bc_std_overall = np.zeros([12,4])
ccn_std_overall = np.zeros([12,4])

i_counter = 0
for counter in ["ss1", "ss2", "ss3", "ss4"]:
    f1 = "data/ccn_std_overall_%s.txt" % counter
    f2 = "data/bc_std_overall_%s.txt" % counter

    ccn_std_overall[:,i_counter] = np.loadtxt(f1) 
    bc_std_overall[:,i_counter] = np.loadtxt(f2) 
    i_counter += 1

x_array = np.arange(-4,2)
help_array = ccn_std_overall[0:6,1]
print help_array.shape, x_array.shape
print help_array

plt.clf()
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
#plt.plot(x_array, ccn_std_overall[0:6,0], 'r-', label = '1K S_c = 0.001%')
plt.plot(x_array, ccn_std_overall[0:6,1], 'g-', label = '1K S_c = 0.01%')
plt.plot(x_array, ccn_std_overall[0:6,2], 'b-', label = '1K S_c = 0.1%')
plt.plot(x_array, ccn_std_overall[0:6,3], 'm-', label = '1K S_c = 0.5%')
#plt.plot(x_array, ccn_std_overall[6:11,0], 'r-', label = '10K S_c = 0.001%')
plt.plot(x_array, ccn_std_overall[6:12,1], 'g--', label = '10K S_c = 0.01%')
plt.plot(x_array, ccn_std_overall[6:12,2], 'b--', label = '10K S_c = 0.1%')
plt.plot(x_array, ccn_std_overall[6:12,3], 'm--', label = '10K S_c = 0.5%')

plt.legend(loc = 'upper right')
plt.xlabel("weighting exponent ")
plt.ylabel("averaged std")
fig = plt.gcf()
fig.savefig("figs/total_std_ccn.pdf")

plt.clf()
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
#plt.plot(x_array, ccn_std_overall[0:6,0], 'r-', label = '1K S_c = 0.001%')
plt.plot(x_array, bc_std_overall[0:6,1], 'g-', label = '1K S_c = 0.01%')
plt.plot(x_array, bc_std_overall[0:6,2], 'b-', label = '1K S_c = 0.1%')
plt.plot(x_array, bc_std_overall[0:6,3], 'm-', label = '1K S_c = 0.5%')
#plt.plot(x_array, ccn_std_overall[6:11,0], 'r-', label = '10K S_c = 0.001%')
plt.plot(x_array, bc_std_overall[6:12,1], 'g--', label = '10K S_c = 0.01%')
plt.plot(x_array, bc_std_overall[6:12,2], 'b--', label = '10K S_c = 0.1%')
plt.plot(x_array, bc_std_overall[6:12,3], 'm--', label = '10K S_c = 0.5%')

plt.legend(loc = 'upper right')
plt.xlabel("weighting exponent ")
plt.ylabel("averaged std")
fig = plt.gcf()
fig.savefig("figs/total_std_bc.pdf")



    



