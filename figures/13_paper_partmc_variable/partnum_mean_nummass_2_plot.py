#!/usr/bin/env python

import os, sys
import config
import scipy.io
import scipy.stats
import numpy as np

sys.path.append("../../tool")
import partmc
import mpl_helper
import matplotlib

#calculation of confidence interval

c = config.c_value
n = config.i_loop_max
r = scipy.stats.t.ppf((1 + c) / 2, n - 1)
conf_factor = r / np.sqrt(n)

num_avg_overall = np.zeros([21,config.i_loop_max])
mass_avg_overall = np.zeros([21,config.i_loop_max])
num_std_overall = np.zeros([21,config.i_loop_max])
mass_std_overall = np.zeros([21,config.i_loop_max])

error_ratio_1K_num = np.zeros([7])
error_ratio_1K_mass = np.zeros([7])

error_ratio_10K_num = np.zeros([7])
error_ratio_10K_mass = np.zeros([7])

i_counter = 0
for counter in ["1K_wei\\+1","10K_wei\\+1", "100K_wei\\+1",
                "1K_flat", "10K_flat", "100K_flat",
                "1K_wei-1", "10K_wei-1", "100K_wei-1",
                "1K_wei-2", "10K_wei-2", "100K_wei-2",
                "1K_wei-3", "10K_wei-3", "100K_wei-3",
                "1K_wei-4", "10K_wei-4", "100K_wei-4",
                "1K_mfa", "10K_mfa", "100K_mfa", ]:
    f1 = "data/ensemble_number_%s.txt" % counter
    f2 = "data/ensemble_mass_%s.txt" % counter
    f3 = "data/ensemble_number_std_%s.txt" % counter
    f4 = "data/ensemble_mass_std_%s.txt" % counter

    num_avg_overall[i_counter,:] = np.loadtxt(f1) / 1e6
    mass_avg_overall[i_counter,:] = np.loadtxt(f2)
    num_std_overall[i_counter,:] = np.loadtxt(f3) / 1e6
    mass_std_overall[i_counter,:] = np.loadtxt(f4)

    i_counter +=1

x_array = [1000, 10000, 1e5]
i_weight = 1

(figure, axes_array) = mpl_helper.make_fig_array(2,1, figure_width=config.figure_width_single, 
                                                 left_margin=0.75, vert_sep=0.2)

axes = axes_array[1][0]

axes.errorbar(x_array[0], num_avg_overall[i_weight*3,99], r*num_std_overall[i_weight*3,99], ecolor='g')
axes.errorbar(x_array[0], num_avg_overall[i_weight*3,99], conf_factor * num_std_overall[i_weight*3,99],
                 marker='_', mec = 'k', ms=7, elinewidth = 7, capsize = 0, ecolor='r')
axes.errorbar(x_array[1], num_avg_overall[1+i_weight*3,99], r*num_std_overall[1+i_weight*3,99], ecolor='g')
axes.errorbar(x_array[1], num_avg_overall[1+i_weight*3,99], conf_factor * num_std_overall[1+i_weight*3,99],
                 marker='_', mec = 'k', ms=7, elinewidth = 7, capsize = 0, ecolor='r')
axes.errorbar(x_array[2], num_avg_overall[2+i_weight*3,99], r*num_std_overall[2+i_weight*3,99], ecolor='g')
axes.errorbar(x_array[2], num_avg_overall[2+i_weight*3,99], conf_factor * num_std_overall[2+i_weight*3,99],
                 marker='_', mec = 'k', ms=7, elinewidth = 7, capsize = 0, ecolor='r')
axes.set_xscale("log")
axes.set_yscale("linear")

axes.set_xlim([500, 2e5])
axes.grid(True)

axes.set_ylabel(r"$\overline{\langle N(t) \rangle}$ / $\rm cm^{-3}$")

axes = axes_array[0][0]
axes.errorbar(x_array[0], mass_avg_overall[i_weight*3,99], r*mass_std_overall[i_weight*3,99],  ecolor='g')
axes.errorbar(x_array[0], mass_avg_overall[i_weight*3,99], conf_factor * mass_std_overall[i_weight*3,99],
                  marker='_', mec = 'k', ms= 7, elinewidth = 7, capsize = 0, ecolor='r')
axes.errorbar(x_array[1], mass_avg_overall[1+i_weight*3,99], r*mass_std_overall[1+i_weight*3,99], ecolor='g')
axes.errorbar(x_array[1], mass_avg_overall[1+i_weight*3,99], conf_factor * mass_std_overall[1+i_weight*3,99],
                 marker='_', mec = 'k', ms= 7, elinewidth = 7, capsize = 0, ecolor='r')
axes.errorbar(x_array[2], mass_avg_overall[2+i_weight*3,99], r*mass_std_overall[2+i_weight*3,99], ecolor='g')
axes.errorbar(x_array[2], mass_avg_overall[2+i_weight*3,99], conf_factor * mass_std_overall[2+i_weight*3,99],
                 marker='_', mec = 'k', ms= 7, elinewidth = 7, capsize = 0, ecolor='r')
axes.set_xscale("log")
axes.set_yscale("linear")

axes.set_xlim([500, 2e5])
axes.grid(True)

axes.set_xlabel(r"particle number $N_{\rm p}$")
axes.set_ylabel(r"$\overline{\langle M(t) \rangle}$ / $\rm \mu g \, m^{-3}$")

mpl_helper.remove_fig_array_axes(axes_array)

figure.savefig("figs/partnum_mean_massnum.pdf")

