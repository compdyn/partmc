#!/usr/bin/env python

import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt

for name in ["32_dist_dist",
             "32_local1_dist",
             "32_local2_dist",
             "32_local3_dist",
             "32_none_dist",
             ]:
    filename = "times_%s.txt" % name
    data = np.loadtxt(filename)
    plt.clf()
    plt.plot(data[:,0], data[:,1], label=("%s cpu" % name))
    plt.hold(True)
    plt.plot(data[:,0], data[:,2], label=("%s wall" % name))
    plt.hold(True)
    plt.xlabel("cores")
    plt.ylabel("time (s)")
    plt.xlim([0, 1024])
    plt.legend(loc="lower right")
    fig = plt.gcf()
    fig.savefig("%s.pdf" % name)
