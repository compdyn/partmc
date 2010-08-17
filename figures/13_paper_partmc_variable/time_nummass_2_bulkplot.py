#!/usr/bin/env python

import os, sys
import config
import config_matplotlib
import scipy.io
import numpy as np
import matplotlib

import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc

time_array = np.loadtxt("data/time_array_1K_flat.txt")
array_num = np.loadtxt("data/array_number_1K_flat.txt")
num_avg = np.loadtxt("data/average_number_1K_flat.txt")
num_std = np.loadtxt("data/std_number_1K_flat.txt")
array_mass = np.loadtxt("data/array_mass_1K_flat.txt")
mass_avg = np.loadtxt("data/average_mass_1K_flat.txt")
mass_std = np.loadtxt("data/std_mass_1K_flat.txt")

num_avg = num_avg / 1e6
num_std = num_std / 1e6

print "time, num, std ", time_array, num_avg, num_std

(figure, axes, cbar_axes) = config_matplotlib.make_fig(kind="1d")
axes.errorbar(time_array, num_avg, num_std, fmt='b-')
axes.set_ylabel(r"number conc. / $\rm cm^{-3}$")
axes.set_ybound(lower=0)

axes.grid(True)
#axes.grid(True, which='minor')
axes.minorticks_on()

axes2 = axes.twinx()
axes2.errorbar(time_array, mass_avg, mass_std, fmt='r-')
axes2.set_ylabel(r"mass conc. / $(\rm \mu g\ m^{-3})$")
axes2.set_ybound(lower=0)

axes.set_xlim(0, 24)
axes.set_xticks([0, 6, 12, 18, 24])
axes.set_xlabel(r"time elapsed / h")

config_matplotlib.label_plot_line(axes, time_array, num_avg, 6, "number",
                                  verticalalignment="bottom",
                                  horizontalalignment="right")
config_matplotlib.label_plot_line(axes2, time_array, mass_avg, 8, "mass",
                                  verticalalignment="top",
                                  horizontalalignment="left")

figure.savefig("figs/time_nummass_bulk.pdf")


(figure, axes, cbar_axes) = config_matplotlib.make_fig(kind="1d")
axes.errorbar(time_array, num_std/num_avg, fmt='b-')
axes.set_ylabel(r"${\rm CV}(N(t))$")
axes.set_ybound(lower=0)

axes.grid(True)
#axes.grid(True, which='minor')
axes.minorticks_on()

axes2 = axes.twinx()
axes2.errorbar(time_array, mass_std/mass_avg, fmt='r-')
axes2.set_ylabel(r"${\rm CV}(M(t))$")
axes2.set_ybound(lower=0)

axes.set_xlim(0, 24)
axes.set_xticks([0, 6, 12, 18, 24])
axes.set_xlabel(r"time elapsed / h")

config_matplotlib.label_plot_line(axes, time_array, num_std/num_avg, 14, "number",
                                  verticalalignment="bottom",
                                  horizontalalignment="right")
config_matplotlib.label_plot_line(axes2, time_array, mass_std/mass_avg, 12, "mass",
                                  verticalalignment="top",
                                  horizontalalignment="left")

figure.savefig("figs/time_cv_nummass_bulk.pdf")

