#!/usr/bin/env python

import os, sys
import scipy.io
import numpy as np

sys.path.append("../../tool")
import partmc
import mpl_helper
import matplotlib
import config

time_array = np.loadtxt("data/time_array_1K_flat.txt")
array_num = np.loadtxt("data/array_number_1K_flat.txt")
num_avg = np.loadtxt("data/average_number_1K_flat.txt")
num_std = np.loadtxt("data/std_number_1K_flat.txt")
array_mass = np.loadtxt("data/array_mass_1K_flat.txt")
mass_avg = np.loadtxt("data/average_mass_1K_flat.txt")
mass_std = np.loadtxt("data/std_mass_1K_flat.txt")

num_avg = num_avg / 1e6
num_std = num_std / 1e6

(figure, axes_array) = mpl_helper.make_fig_array(1,2, figure_width=config.figure_width_double, 
                                                 left_margin=0.7, right_margin=0.6, horiz_sep=1,share_y_axes=False)

axes = axes_array[0][0]

axes.errorbar(time_array, num_avg, num_std, fmt='b-')
axes.set_ylabel(r"total number $N(t)$ / $\rm cm^{-3}$")
axes.set_ybound(lower=0)
axes.set_yticks([0, 4000, 8000, 12000, 16000])
axes.grid(True)

axes2 = axes.twinx()
axes2.errorbar(time_array, mass_avg, mass_std, fmt='r-')
axes2.set_ylabel(r"total mass $M(t)$ / $(\rm \upmu g\ m^{-3})$")
axes2.set_ybound(lower=0)
axes2.set_yticks([0, 10, 20, 30, 40])

axes.set_xlim(0, 24)
axes.set_xticks([0, 6, 12, 18, 24])
axes.set_xlabel(r"time elapsed / h")

mpl_helper.label_plot_line(axes, time_array, num_avg, 8, "number",
                                  verticalalignment="bottom",
                                  horizontalalignment="right")
mpl_helper.label_plot_line(axes2, time_array, mass_avg, 8, "mass",
                                  verticalalignment="top",
                                  horizontalalignment="left")

axes = axes_array[0][1]
axes.errorbar(time_array, num_std/num_avg, fmt='b-')
axes.set_ylabel(r"${\rm CV}(N(t))$")
axes.set_ylim(0, 0.04)
axes.set_yticks([0, 0.01, 0.02, 0.03, 0.04])
axes.grid(True)

axes2 = axes.twinx()
axes2.errorbar(time_array, mass_std/mass_avg, fmt='r-')
axes2.set_ylabel(r"${\rm CV}(M(t))$")
axes2.set_ylim(0, 0.2)
axes2.set_yticks([0, 0.05, 0.1, 0.15, 0.2])

axes.set_xlim(0, 24)
axes.set_xticks([0, 6, 12, 18, 24])
axes.set_xlabel(r"time elapsed / h")

mpl_helper.label_plot_line(axes, time_array, num_std/num_avg, 14, "number",
                                  verticalalignment="bottom",
                                  horizontalalignment="right")
mpl_helper.label_plot_line(axes2, time_array, mass_std/mass_avg, 12, "mass",
                                  verticalalignment="top",
                                  horizontalalignment="left")

figure.savefig("figs/time_cv_nummass_bulk.pdf")

