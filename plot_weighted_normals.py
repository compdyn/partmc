#!/usr/bin/env python

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt

x = np.linspace(-8, 2, 100)
N = 100
s = 1.5
xb = -5
al = -1
x0 = -4

n = N / np.sqrt(2 * np.pi * s**2) * np.exp(-(x - xb)**2 / (2*s**2))
w = 10**(al*(x-x0))
c0 = n / w

xbp = xb - np.log(10) * al * s**2
Np = N * 10**(al*x0) * np.exp((xbp**2 - xb**2) / (2*s**2))

c1 = Np / np.sqrt(2 * np.pi * s**2) * np.exp(-(x - xbp)**2 / (2*s**2))

plt.semilogy(x,n)
plt.semilogy(x,w)
plt.semilogy(x,c0)
plt.semilogy(x,c1,'x')
plt.grid(True)
plt.legend(('n', 'w', 'c0', 'c1'))
plt.savefig("weighted_normals.pdf")
