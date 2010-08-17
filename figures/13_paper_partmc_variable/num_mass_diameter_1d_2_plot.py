#!/usr/bin/env python2.5

import os, sys
import config
import scipy.io
import numpy as np

sys.path.append("../../tool")
import partmc
import mpl_helper
import matplotlib

num_avg = np.zeros([3])
mass_avg = np.zeros([3])
num_std = np.zeros([3])
mass_std = np.zeros([3])

x_array = np.loadtxt("data/1d_10K_wei+1_12_x_values.txt")
num_avg1 = np.loadtxt("data/1d_10K_wei+1_12_hist_array_gav_num.txt")
num_avg2 = np.loadtxt("data/1d_10K_wei-1_12_hist_array_gav_num.txt")
num_avg3 = np.loadtxt("data/1d_10K_wei-4_12_hist_array_gav_num.txt")

num_std1 = np.loadtxt("data/1d_10K_wei+1_12_e_bars_num.txt")
num_std2 = np.loadtxt("data/1d_10K_wei-1_12_e_bars_num.txt")
num_std3 = np.loadtxt("data/1d_10K_wei-4_12_e_bars_num.txt")

mass_avg1 = np.loadtxt("data/1d_10K_wei+1_12_hist_array_gav_mass.txt")
mass_avg2 = np.loadtxt("data/1d_10K_wei-1_12_hist_array_gav_mass.txt")
mass_avg3 = np.loadtxt("data/1d_10K_wei-4_12_hist_array_gav_mass.txt")

mass_std1 = np.loadtxt("data/1d_10K_wei+1_12_e_bars_mass.txt")
mass_std2 = np.loadtxt("data/1d_10K_wei-1_12_e_bars_mass.txt")
mass_std3 = np.loadtxt("data/1d_10K_wei-4_12_e_bars_mass.txt")

x_array = x_array * 1e6
num_avg1 = num_avg1 / 1e6
num_std1 = num_std1 / 1e6
mass_avg1 = mass_avg1 * 1e9
mass_std1 = mass_std1 * 1e9

(figure, axes_array) = mpl_helper.make_fig_array(3,2, figure_width=config.figure_width_double, 
                                                 left_margin=0.7, right_margin=0.6, vert_sep=0.3)

axes = axes_array[2][0]
print 'num1 ', x_array, num_avg1, num_std1
axes.errorbar(x_array, num_avg1, num_std1, fmt='g-')
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_ylabel(r"$n(D)$ / $\rm cm^{-3}$")
axes.set_ylim(1e-2, 1e5)
axes.set_xlim(5e-3, 5)
axes.grid(True)

axes = axes_array[1][0]
axes.errorbar(x_array, num_avg2, num_std2, fmt='g-')
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_ylabel(r"$n(D)$ / $\rm cm^{-3}$")
#axes.set_ylim(1e-2, 1e5)
axes.set_xlim(5e-3, 5)
axes.grid(True)

axes = axes_array[0][0]
axes.errorbar(x_array, num_avg3, num_std3, fmt='g-')
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_ylabel(r"$n(D)$ / $\rm cm^{-3}$")
#axes.set_ylim(1e-2, 1e5)
axes.set_xlim(5e-3, 5)
axes.grid(True)
axes.set_xlabel(r"diameter / $\rm \mu m$")

axes = axes_array[2][1]
axes.errorbar(x_array, mass_avg1,mass_std1, fmt='b-')
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_ylabel(r"$m(D)$ / $\rm \mu g \, m^{-3}$")
#axes.set_xlim(5e-3, 5)
axes.set_ylim(1e-4, 1e2)
axes.grid(True)

axes = axes_array[1][1]
axes.errorbar(x_array, mass_avg2, mass_std2, fmt='b-')
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_ylabel(r"$m(D)$ / $\rm \mu g \, m^{-3}$")
#axes.set_ylim(1e-4, 1e2)
axes.set_xlim(5e-3, 5)
axes.grid(True)

axes = axes_array[0][1]
axes.errorbar(x_array, mass_avg3, mass_std3, fmt='b-')
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_ylabel(r"$m(D)$ / $\rm \mu g \, m^{-3}$")
#axes.set_ylim(1e-4, 1e2)
axes.set_xlim(5e-3, 5)
axes.grid(True)
axes.set_xlabel(r"diameter / $\rm \mu m$")

mpl_helper.remove_fig_array_axes(axes_array)

figure.savefig("figs/num_mass_diameter_1d.pdf")

