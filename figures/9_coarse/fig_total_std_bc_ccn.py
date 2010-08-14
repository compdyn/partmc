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

bc_std_overall = np.zeros([12,4]) # 12 for weighting schemes, 4 for ss-values
ccn_std_overall = np.zeros([12,4])

i_counter = 0
for counter in ["ss1", "ss2", "ss3", "ss4"]:
    print "counter ", i_counter, counter
    print "ccn_std_overall", ccn_std_overall
    f1 = "data/ccn_std_overall_%s.txt" % counter
    f2 = "data/bc_std_overall_%s.txt" % counter
    print "f1 ", f1
    print "content ", np.loadtxt(f1)
    ccn_std_overall[:,i_counter] = np.loadtxt(f1) 
    bc_std_overall[:,i_counter] = np.loadtxt(f2) 
    i_counter += 1

x_array = np.arange(1, -5, -1)
help_array1 = ccn_std_overall[0:6,1]
help_array2 = ccn_std_overall[6:12,1]
print help_array1.shape, help_array2.shape, x_array.shape
print ccn_std_overall
print x_array, help_array1, help_array2

plt.clf()
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
plt.plot(x_array, ccn_std_overall[0:6,1], 'g-', label = '1K S_c = 0.01%')
plt.plot(x_array, ccn_std_overall[0:6,2], 'b-', label = '1K S_c = 0.1%')
plt.plot(x_array, ccn_std_overall[0:6,3], 'm-', label = '1K S_c = 0.5%')

#plt.plot(-3, ccn_std_overall[18,1], 'go', label = '1K S_c = 0.01% mfa')
#plt.plot(-3, ccn_std_overall[18,2], 'bo', label = '1K S_c = 0.1% mfa')
#plt.plot(-3, ccn_std_overall[18,3], 'mo', label = '1K S_c = 0.5% mfa')

plt.plot(x_array, ccn_std_overall[6:12,1], 'g--', label = '10K S_c = 0.01%')
plt.plot(x_array, ccn_std_overall[6:12,2], 'b--', label = '10K S_c = 0.1%')
plt.plot(x_array, ccn_std_overall[6:12,3], 'm--', label = '10K S_c = 0.5%')

#plt.plot(-3, ccn_std_overall[19,1], 'go', label = '10K S_c = 0.01% mfa')
#plt.plot(-3, ccn_std_overall[19,2], 'bo', label = '10K S_c = 0.1% mfa')
#plt.plot(-3, ccn_std_overall[19,3], 'mo', label = '10K S_c = 0.5% mfa')

#plt.plot(x_array, ccn_std_overall[12:18,1], 'g--', label = '100K S_c = 0.01%')
#plt.plot(x_array, ccn_std_overall[12:18,2], 'b--', label = '100K S_c = 0.1%')
#plt.plot(x_array, ccn_std_overall[12:18,3], 'm--', label = '100K S_c = 0.5%')

#plt.plot(-3, ccn_std_overall[20,1], 'go', label = '100K S_c = 0.01% mfa')
#plt.plot(-3, ccn_std_overall[20,2], 'bo', label = '100K S_c = 0.1% mfa')
#plt.plot(-3, ccn_std_overall[20,3], 'mo', label = '100K S_c = 0.5% mfa')


plt.grid()
plt.legend(loc = 'upper right')
plt.xlabel("weighting exponent ")
plt.ylabel("averaged std/avg")
fig = plt.gcf()
fig.savefig("figs/total_std_ccn.pdf")

plt.clf()
a = plt.gca()
a.set_xscale("linear")
a.set_yscale("log")
plt.plot(x_array, bc_std_overall[0:6,1], 'g-', label = '1K S_c = 0.01%')
plt.plot(x_array, bc_std_overall[0:6,2], 'b-', label = '1K S_c = 0.1%')
plt.plot(x_array, bc_std_overall[0:6,3], 'm-', label = '1K S_c = 0.5%')

#plt.plot(-3, bc_std_overall[18,1], 'go', label = '1K S_c = 0.01% mfa')
#plt.plot(-3, bc_std_overall[18,2], 'bo', label = '1K S_c = 0.1% mfa')
#plt.plot(-3,bcn_std_overall[18,3], 'mo', label = '1K S_c = 0.5% mfa')

plt.plot(x_array, bc_std_overall[6:12,1], 'g--', label = '10K S_c = 0.01%')
plt.plot(x_array, bc_std_overall[6:12,2], 'b--', label = '10K S_c = 0.1%')
plt.plot(x_array, bc_std_overall[6:12,3], 'm--', label = '10K S_c = 0.5%')

#plt.plot(-3, bc_std_overall[19,1], 'go', label = '10K S_c = 0.01% mfa')
#plt.plot(-3, bc_std_overall[19,2], 'bo', label = '10K S_c = 0.1% mfa')
#plt.plot(-3, bc_std_overall[19,3], 'mo', label = '10K S_c = 0.5% mfa')

#plt.plot(x_array, bc_std_overall[12:18,1], 'g--', label = '100K S_c = 0.01%')
#plt.plot(x_array, bc_std_overall[12:18,2], 'b--', label = '100K S_c = 0.1%')
#plt.plot(x_array, bc_std_overall[12:18,3], 'm--', label = '100K S_c = 0.5%')

#plt.plot(-3, bc_std_overall[20,1], 'go', label = '100K S_c = 0.01% mfa')
#plt.plot(-3, bc_std_overall[20,2], 'bo', label = '100K S_c = 0.1% mfa')
#plt.plot(-3, bc_std_overall[20,3], 'mo', label = '100K S_c = 0.5% mfa')

plt.grid()
plt.legend(loc = 'upper right')
plt.xlabel("weighting exponent ")
plt.ylabel("averaged std/avg")
fig = plt.gcf()
fig.savefig("figs/total_std_bc.pdf")

plt.clf()
a = plt.gca()
a.set_xscale("log")
a.set_yscale("log")
plt.plot(bc_std_overall[0:6,1], ccn_std_overall[0:6,1],'g-x', label = '1K S_c = 0.01%')
plt.plot(bc_std_overall[0:6,2], ccn_std_overall[0:6,2], 'b-x', label = '1K S_c = 0.1%')
plt.plot(bc_std_overall[0:6,3], ccn_std_overall[0:6,3], 'm-x', label = '1K S_c = 0.5%')

#plt.plot(bc_std_overall[18,1], ccn_std_overall[18,1], 'go', label = '1K S_c = 0.01% mfa')
#plt.plot(bc_std_overall[18,2], ccn_std_overall[18,2], 'bo', label = '1K S_c = 0.1% mfa')
#plt.plot(bc_std_overall[18,3], ccn_std_overall[18,3], 'mo', label = '1K S_c = 0.5% mfa')

plt.plot(bc_std_overall[6:12,1], ccn_std_overall[6:12,1], 'g--x', label = '10K S_c = 0.01%')
plt.plot(bc_std_overall[6:12,2], ccn_std_overall[6:12,2], 'b--x', label = '10K S_c = 0.1%')
plt.plot(bc_std_overall[6:12,3], ccn_std_overall[6:12,3], 'm--x', label = '10K S_c = 0.5%')

#plt.plot(bc_std_overall[19,1], ccn_std_overall[19,1], 'go', label = '10K S_c = 0.01% mfa')
#plt.plot(bc_std_overall[19,2], ccn_std_overall[19,2], 'bo', label = '10K S_c = 0.1% mfa')
#plt.plot(bc_std_overall[19,3], ccn_std_overall[19,3], 'mo', label = '10K S_c = 0.5% mfa')

#plt.plot(bc_std_overall[12:18,1], ccn_std_overall[12:18,1], 'g--', label = '100K S_c = 0.01%')
#plt.plot(bc_std_overall[12:18,2], ccn_std_overall[12:18,2], 'b--', label = '100K S_c = 0.1%')
#plt.plot(bc_std_overall[12:18,3], ccn_std_overall[12:18,3], 'm--', label = '100K S_c = 0.5%')

#plt.plot(bc_std_overall[20,1], ccn_std_overall[20,1], 'go', label = '100K S_c = 0.01% mfa')
#plt.plot(bc_std_overall[20,2], ccn_std_overall[20,2], 'bo', label = '100K S_c = 0.1% mfa')
#plt.plot(bc_std_overall[20,3], ccn_std_overall[20,3], 'mo', label = '100K S_c = 0.5% mfa')

plt.grid()
plt.legend(loc = 'lower right')
plt.xlabel("CV(N) ")
plt.ylabel("CV(M)")
fig = plt.gcf()
fig.savefig("figs/total_std_bc_versus_ccn.pdf")


    



