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

num_avg_overall = np.zeros([config.i_loop_max])
mass_avg_overall = np.zeros([config.i_loop_max])

f1 = "data/ensemble_norm_num_brownian_60s.txt" 
f2 = "data/ensemble_norm_mass_brownian_60s.txt" 
#f3 = "data/ensemble_norm_num_brownian_1200s.txt" 
#f4 = "data/ensemble_norm_mass_brownian_1200s.txt" 

num_avg_overall = np.loadtxt(f1) / 1e6
mass_avg_overall = np.loadtxt(f2)*1e9

#num_avg_overall_6s = np.loadtxt(f3) / 1e6
#mass_avg_overall_6s = np.loadtxt(f4)*1e9

(figure, axes_array) = mpl_helper.make_fig_array(1,2, figure_width=config.figure_width_double, 
                                                 left_margin=0.75, horiz_sep=0.6, share_y_axes = False)

axes = axes_array[0][0]

large_time_step = axes.plot(num_avg_overall, 'b-')
axes.set_xscale("log")
axes.set_yscale("log")

axes.set_xlim([0, 100])
axes.set_ylim([num_avg_overall.min(), num_avg_overall.max()])
axes.grid(True)
axes.set_xlabel(r"ensemble size")
axes.set_ylabel(r" $\| \langle  n_{\rm p}(24\, {\rm h}) - n_{\rm sect}(24 \, {\rm h}) \rangle \|_2$")
#axes.legend((large_time_step, small_time_step), (r"60 s", r"1200 s"))

axes = axes_array[0][1]
large_time_step = axes.plot(mass_avg_overall, 'b-')
axes.set_xscale("log")
axes.set_yscale("log")

axes.set_xlim([0, 100])
axes.set_ylim([mass_avg_overall.min(), mass_avg_overall.max()])
axes.grid(True)

axes.set_xlabel(r"ensemble size")
axes.set_ylabel(r" $\| \langle  m_{\rm p}(24\, {\rm h}) - m_{\rm sect}(24 \, {\rm h}) \rangle \|_2$")
#axes.legend((large_time_step, small_time_step), (r"60 s", r"1200 s"))
figure.savefig("figs/brownian_ensemble_averages.pdf")

