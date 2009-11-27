#!/usr/bin/env python2.5

import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
import config

ccn_cn_array = np.loadtxt("data/ccn_cn.txt")

plt.plot(ccn_cn_array[:,0]/3600., ccn_cn_array[:,1], label = 'S = %.2f' % config.s_crit_1)
plt.plot(ccn_cn_array[:,0]/3600., ccn_cn_array[:,2], label = 'S = %.2f' % config.s_crit_2)
plt.plot(ccn_cn_array[:,0]/3600., ccn_cn_array[:,3], label = 'S = %.2f' % config.s_crit_3)
plt.axis([0, 48, 0, 1])
plt.legend(loc = 'center right')
plt.xlabel("time (hours)")
plt.ylabel("ccn/cn")

fig = plt.gcf()
fig.savefig("figs/ccn_time.pdf")


    



