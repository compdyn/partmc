#!/usr/bin/env python

import os, sys
import scipy.io
import scipy.stats
import numpy as np

sys.path.append("../../tool")
import partmc
import mpl_helper
import matplotlib
import config

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

i_weight = 0
for weight in ["wei+1", "flat", "wei-1", "wei-2", "wei-3", "wei-4", "mfa"]:

    index_1K = 3*i_weight
    index_10K = 3*i_weight + 1

    error_ratio_1K_num[i_weight] = (num_avg_overall[index_1K+1, 99] - num_avg_overall[index_1K, 99])/(r*num_std_overall[index_1K,99])
    error_ratio_10K_num[i_weight] = (num_avg_overall[index_10K+1, 99] - num_avg_overall[index_10K, 99])/(r*num_std_overall[index_10K,99])

    error_ratio_1K_mass[i_weight] = (mass_avg_overall[index_1K+1, 99] - mass_avg_overall[index_1K, 99])/(r*mass_std_overall[index_1K,99])
    error_ratio_10K_mass[i_weight] = (mass_avg_overall[index_10K+1, 99] - mass_avg_overall[index_10K, 99])/(r*mass_std_overall[index_10K,99])
    i_weight += 1

weight_array = [1, 0, -1, -2, -3, -4]

(figure, axes) = mpl_helper.make_fig(figure_width = 0.7*config.figure_width_double, top_margin=0.2, bottom_margin=0.4,left_margin=0.6, right_margin=1.3)
l_n_1e3 = axes.plot(weight_array, error_ratio_1K_num[0:6], label = r"${\rm ER_N}(10^3,10^4)$")
l_n_1e4 = axes.plot(weight_array, error_ratio_10K_num[0:6], label = r"${\rm ER_N}(10^4,10^5)$")
l_m_1e3 = axes.plot(weight_array, error_ratio_1K_mass[0:6], label = r"${\rm ER_M}(10^3,10^4)$")
l_m_1e4 = axes.plot(weight_array, error_ratio_10K_mass[0:6], label = r"${\rm ER_M}(10^4,10^5)$")
axes.grid(True)
#axes.legend()
axes.set_ylim(-0.8, 0.5)
axes.set_xlabel(r"weighting exponent $\alpha$ ")
axes.set_ylabel(r"error ratio ER")
axes.legend((l_n_1e3,l_n_1e4,l_m_1e3,l_m_1e4), 
              (r"${\rm ER_N}(10^3,10^4)$", r"${\rm ER_N}(10^4,10^5)$", 
               r"${\rm ER_M}(10^3,10^4)$", r"${\rm ER_M}(10^4,10^5)$"), loc='center left', bbox_to_anchor = (1.0,0.5),ncol=1)

figure.savefig("figs/er_alpha.pdf")
