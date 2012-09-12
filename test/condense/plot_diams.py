#!/usr/bin/env python

import sys
sys.path.append('../../../../tool')
import mpl_helper
import matplotlib
import numpy

(fig, axes) = mpl_helper.make_fig(5)
data = numpy.loadtxt("diams.txt")
for i in range(1, data.shape[1]):
    axes.semilogy(data[:,0], data[:,i])

fig.savefig("diams.pdf")
