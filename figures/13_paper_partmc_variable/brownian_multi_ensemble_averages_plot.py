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

c = config.c_value
n = config.i_ens_max
r = scipy.stats.t.ppf((1 + c) / 2, n - 1)
conf_factor = r / np.sqrt(n)


f1 = "data/avg_ensemble_norm_num_brownian_600s.txt"
f2 = "data/avg_ensemble_norm_mass_brownian_600s.txt" 
f3 = "data/std_ensemble_norm_num_brownian_600s.txt"
f4 = "data/std_ensemble_norm_mass_brownian_600s.txt" 

f5 = "data/avg_ensemble_norm_num_brownian_60s.txt"
f6 = "data/avg_ensemble_norm_mass_brownian_60s.txt" 
f7 = "data/std_ensemble_norm_num_brownian_60s.txt"
f8 = "data/std_ensemble_norm_mass_brownian_60s.txt" 

f9 = "data/avg_ensemble_norm_num_brownian_6s.txt"
f10 = "data/avg_ensemble_norm_mass_brownian_6s.txt" 
f11 = "data/std_ensemble_norm_num_brownian_6s.txt"
f12 = "data/std_ensemble_norm_mass_brownian_6s.txt" 

f13 = "data/avg_ensemble_norm_num_brownian_1200s.txt"
f14 = "data/avg_ensemble_norm_mass_brownian_1200s.txt" 
f15 = "data/std_ensemble_norm_num_brownian_1200s.txt"
f16 = "data/std_ensemble_norm_mass_brownian_1200s.txt" 


num_avg_overall_600s = np.loadtxt(f1) / 1e6
mass_avg_overall_600s = np.loadtxt(f2)*1e9
num_std_overall_600s = np.loadtxt(f3) / 1e6
mass_std_overall_600s = np.loadtxt(f4)*1e9

num_avg_overall_60s = np.loadtxt(f5) / 1e6
mass_avg_overall_60s = np.loadtxt(f6)*1e9
num_std_overall_60s = np.loadtxt(f7) / 1e6
mass_std_overall_60s = np.loadtxt(f8)*1e9

num_avg_overall_6s = np.loadtxt(f9) / 1e6
mass_avg_overall_6s = np.loadtxt(f10)*1e9
num_std_overall_6s = np.loadtxt(f11) / 1e6
mass_std_overall_6s = np.loadtxt(f12)*1e9

num_avg_overall_1200s = np.loadtxt(f13) / 1e6
mass_avg_overall_1200s = np.loadtxt(f14)*1e9
num_std_overall_1200s = np.loadtxt(f15) / 1e6
mass_std_overall_1200s = np.loadtxt(f16)*1e9

x_array = range(1,101)

(figure, axes_array) = mpl_helper.make_fig_array(1,2, figure_width=config.figure_width_double, 
                                                 left_margin=0.75, horiz_sep=0.6, share_y_axes = False)

axes = axes_array[0][0]

#time_step_6s = axes.errorbar(x_array, num_avg_overall_6s, conf_factor*num_std_overall_6s, fmt = 'g-')
#time_step_60s = axes.errorbar(x_array, num_avg_overall_60s, conf_factor*num_std_overall_60s, fmt = 'b-')
#time_step_600s = axes.errorbar(x_array, num_avg_overall_600s, conf_factor*num_std_overall_600s, fmt = 'r-')
#time_step_1200s = axes.errorbar(x_array, num_avg_overall_1200s, conf_factor*num_std_overall_1200s, fmt = 'm-')

time_step_6s = axes.plot(x_array, num_avg_overall_6s,  'g-')
time_step_6s_high = axes.plot(x_array, num_avg_overall_6s+conf_factor*num_std_overall_6s,  'g:', linewidth=0.2)
time_step_6s_low = axes.plot(x_array, num_avg_overall_6s-conf_factor*num_std_overall_6s,  'g:', linewidth=0.2)
time_step_60s = axes.plot(x_array, num_avg_overall_60s, 'b-')
time_step_60s_high = axes.plot(x_array, num_avg_overall_60s+conf_factor*num_std_overall_60s,  'b:', linewidth=0.2)
time_step_60s_low = axes.plot(x_array, num_avg_overall_60s-conf_factor*num_std_overall_60s,  'b:', linewidth=0.2)
time_step_600s = axes.plot(x_array, num_avg_overall_600s, 'r-')
time_step_600s_high = axes.plot(x_array, num_avg_overall_600s+conf_factor*num_std_overall_600s,  'r:', linewidth=0.2)
time_step_600s_low = axes.plot(x_array, num_avg_overall_600s-conf_factor*num_std_overall_600s,  'r:', linewidth=0.2)
time_step_1200s = axes.plot(x_array, num_avg_overall_1200s, 'm-')
time_step_1200s_high = axes.plot(x_array, num_avg_overall_1200s+conf_factor*num_std_overall_1200s,  'm:', linewidth=0.2)
time_step_1200s_low = axes.plot(x_array, num_avg_overall_1200s-conf_factor*num_std_overall_1200s,  'm:', linewidth=0.2)
axes.set_xscale("log")
axes.set_yscale("log")

axes.set_xlim([0, 100])
axes.set_ylim([num_avg_overall_600s.min(), num_avg_overall_600s.max()])
axes.grid(True)
axes.set_xlabel(r"ensemble size")
axes.set_ylabel(r" $\| \langle  n_{\rm p}(24\, {\rm h}) - n_{\rm sect}(24 \, {\rm h}) \rangle \|_2$")
axes.legend((time_step_6s, time_step_60s, time_step_600s, time_step_1200s), (r"6 s", r"60 s", r"600 s", r"1200 s"))

axes = axes_array[0][1]

time_step_6s = axes.plot(x_array, mass_avg_overall_6s,  'g-')
time_step_6s_high = axes.plot(x_array, mass_avg_overall_6s+conf_factor*mass_std_overall_6s,  'g:', linewidth=0.2)
time_step_6s_low = axes.plot(x_array, mass_avg_overall_6s-conf_factor*mass_std_overall_6s,  'g:', linewidth=0.2)
time_step_60s = axes.plot(x_array, mass_avg_overall_60s, 'b-')
time_step_60s_high = axes.plot(x_array, mass_avg_overall_60s+conf_factor*mass_std_overall_60s,  'b:', linewidth=0.2)
time_step_60s_low = axes.plot(x_array, mass_avg_overall_60s-conf_factor*mass_std_overall_60s,  'b:', linewidth=0.2)
time_step_600s = axes.plot(x_array, mass_avg_overall_600s, 'r-')
time_step_600s_high = axes.plot(x_array, mass_avg_overall_600s+conf_factor*mass_std_overall_600s,  'r:', linewidth=0.2)
time_step_600s_low = axes.plot(x_array, mass_avg_overall_600s-conf_factor*mass_std_overall_600s,  'r:', linewidth=0.2)
time_step_1200s = axes.plot(x_array, mass_avg_overall_1200s, 'm-')
time_step_1200s_high = axes.plot(x_array, mass_avg_overall_1200s+conf_factor*mass_std_overall_1200s,  'm:', linewidth=0.2)
time_step_1200s_low = axes.plot(x_array, mass_avg_overall_1200s-conf_factor*mass_std_overall_1200s,  'm:', linewidth=0.2)
axes.set_xscale("log")
axes.set_yscale("log")

axes.set_xlim([0, 100])
axes.set_ylim([mass_avg_overall_600s.min(), mass_avg_overall_600s.max()])
axes.grid(True)

axes.set_xlabel(r"ensemble size")
axes.set_ylabel(r" $\| \langle  m_{\rm p}(24\, {\rm h}) - m_{\rm sect}(24 \, {\rm h}) \rangle \|_2$")
axes.legend((time_step_6s, time_step_60s, time_step_600s, time_step_1200s), (r"6 s", r"60 s", r"600 s", r"1200 s"))

figure.savefig("figs/brownian_multi_ensemble_averages.pdf")

