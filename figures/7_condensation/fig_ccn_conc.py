#!/usr/bin/env python2.5

import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
import config

ccn_array_wc = np.loadtxt("data/ccn_conc_wc.txt")
ccn_array_nc = np.loadtxt("data/ccn_conc_nc.txt")

plt.plot(ccn_array_wc[:,0]/3600., ccn_array_wc[:,1], 'b-', label = 'S = %.2f wc' % config.s_crit_1)
plt.plot(ccn_array_wc[:,0]/3600., ccn_array_wc[:,2], 'g-', label = 'S = %.2f wc' % config.s_crit_2)
plt.plot(ccn_array_wc[:,0]/3600., ccn_array_wc[:,3], 'r-', label = 'S = %.2f wc' % config.s_crit_3)

plt.plot(ccn_array_nc[:,0]/3600., ccn_array_nc[:,1], 'b--', label = 'S = %.2f nc' % config.s_crit_1)
plt.plot(ccn_array_nc[:,0]/3600., ccn_array_nc[:,2], 'g--', label = 'S = %.2f nc' % config.s_crit_2)
plt.plot(ccn_array_nc[:,0]/3600., ccn_array_nc[:,3], 'r--', label = 'S = %.2f nc' % config.s_crit_3)

#plt.axis([0, 24, 0, 1])
plt.legend(loc = 'upper right')
plt.xlabel("time (hours)")
plt.ylabel("ccn concentration (m^{-3})")

fig = plt.gcf()
fig.savefig("figs/ccn_conc.pdf")


    



