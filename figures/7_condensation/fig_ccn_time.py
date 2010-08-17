#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
import config

ccn_cn_array = np.loadtxt("data/ccn_cn_laura_100_numbers.txt")
plt.plot(ccn_cn_array[:,0]/3600., ccn_cn_array[:,4], label = 'all')
plt.plot(ccn_cn_array[:,0]/3600., ccn_cn_array[:,1], label = 'S = %.2f' % config.s_crit_1)
plt.plot(ccn_cn_array[:,0]/3600., ccn_cn_array[:,2], label = 'S = %.2f' % config.s_crit_2)
plt.plot(ccn_cn_array[:,0]/3600., ccn_cn_array[:,3], label = 'S = %.2f' % config.s_crit_3)
#plt.axis([0, 24, 0, 1])
plt.legend(loc = 'center right')
plt.title("data/scaleBC3/rate_100")
plt.grid(True)
plt.xlabel("time (hours)")
plt.ylabel("ccn in m^{-3}")

fig = plt.gcf()
fig.savefig("figs/ccn_time_laura_100_numbers.pdf")


    



