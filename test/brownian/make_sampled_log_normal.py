#!/usr/bin/env python

import math

D_min = 5e-9
D_max = 5e-7
n_bin = 20

num_conc = 1e11
geom_mean_diam = 5e-8
log10_geom_std_dev = 0.24

def cdf(x):
    return 0.5 * (1 + math.erf((x - math.log10(geom_mean_diam)) \
                                   / (math.sqrt(2) * log10_geom_std_dev)))

edges = []
num_concs = []

for i_bin in range(n_bin):
    x_min = math.log10(D_min)
    x_max = math.log10(D_max)
    x_0 = float(i_bin) / n_bin * (x_max - x_min) + x_min
    x_1 = float(i_bin + 1) / n_bin * (x_max - x_min) + x_min
    if (i_bin == 0):
        edges.append(10**(x_0))
    edges.append(10**(x_1))
    num_concs.append(num_conc * (cdf(x_1) - cdf(x_0)))

f = open("aero_init_size_dist.dat", "w")
f.write("# first line is diam of bin edges (m)\n")
f.write("# second line is num_conc in each bin (m^{-3})\n")
f.write("diam %s\n" % " ".join([str(x) for x in edges]))
f.write("num_conc %s\n" % " ".join([str(x) for x in num_concs]))
f.close()
