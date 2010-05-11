#!/usr/bin/env python

import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt

for name in ["central_dist",
             "collect_dist",
             "dist_central",
             "dist_dist",
             "dist_none",
             "dist_single",
             "local1_dist",
             "local2_dist",
             "local3_dist",
             ]:
    filename = "times_%s.txt" % name
    data = np.loadtxt(filename)
    plt.clf()
    plt.plot(data[:,0] * 8, data[:,1], label=("%s cpu" % name))
    plt.hold(True)
    plt.plot(data[:,0] * 8, data[:,2], label=("%s wall" % name))
    plt.hold(True)
    plt.xlabel("cores")
    plt.ylabel("time (s)")
    plt.xlim([0, 128])
    plt.legend(loc="lower right")
    fig = plt.gcf()
    fig.savefig("%s.pdf" % name)
