#!/usr/bin/env python

import os, sys, math
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import matplotlib.collections as mcoll
import numpy as np

figure = plt.figure()
axes = figure.add_subplot(111)
verts = np.array([[1,1], [1,2], [2,2], [2,1]])
c = mcoll.PolyCollection([verts], linewidths = 40)
axes.add_collection(c)
axes.set_xbound(0, 3)
axes.set_ybound(0, 3)
figure.savefig("closed_poly.pdf")
