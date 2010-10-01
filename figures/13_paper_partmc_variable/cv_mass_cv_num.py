#!/usr/bin/env python

import os, sys
import scipy.io
import numpy as np

sys.path.append("../../tool")
import partmc
import mpl_helper
import matplotlib
import config

i_counter = 0
i_counter_p = 0
array_num_std_avg = np.zeros([3,7]) # 3 for particle number cases, 7 for weighting cases
array_mass_std_avg = np.zeros([3,7])

num_avg = np.zeros([3,7,25])
num_std = np.zeros([3,7,25])
mass_avg = np.zeros([3,7,25])
mass_std = np.zeros([3,7,25])

for counter_p in ["1K", "10K", "100K"]:
    print counter_p, i_counter_p
    i_counter = 0
    for counter in ["wei\+1", "flat", "wei-1", "wei-2", "wei-3", "wei-4", "mfa"]:
        print counter, i_counter
        f1 = "data/average_number_%s_%s.txt" % (counter_p, counter)
        f2 = "data/std_number_%s_%s.txt" % (counter_p, counter)
        f3 = "data/average_mass_%s_%s.txt" % (counter_p, counter)
        f4 = "data/std_mass_%s_%s.txt" % (counter_p, counter)

        print f1, f2, f3, f4
        num_avg[i_counter_p,i_counter,:] = np.loadtxt(f1)
        num_std[i_counter_p,i_counter,:] = np.loadtxt(f2)
        mass_avg[i_counter_p,i_counter,:] = np.loadtxt(f3)
        mass_std[i_counter_p,i_counter,:] = np.loadtxt(f4)

# average std over time

        array_num_std_avg[i_counter_p,i_counter] = np.average(num_std[i_counter_p,i_counter,:] / 
                                                              num_avg[i_counter_p,i_counter,:])
        array_mass_std_avg[i_counter_p,i_counter] = np.average(mass_std[i_counter_p, i_counter,:] / 
                                                               mass_avg[i_counter_p, i_counter,:])
        i_counter = i_counter + 1
    i_counter_p = i_counter_p + 1

x_array = [1, 0, -1, -2, -3, -4]

(figure, axes) = mpl_helper.make_fig(figure_width=config.figure_width_single, axis_ratio=1)

axes.set_xscale("log")
axes.set_yscale("log")
axes.plot(array_num_std_avg[0,0:6], array_mass_std_avg[0,0:6], 'xr-')
axes.plot(array_num_std_avg[0,6], array_mass_std_avg[0,6], 'or', label = '1K mfa')
axes.plot(array_num_std_avg[1,0:6], array_mass_std_avg[1,0:6], 'xg-')
axes.plot(array_num_std_avg[1,6], array_mass_std_avg[1,6], 'og', label = '10K mfa')
axes.plot(array_num_std_avg[2,0:6], array_mass_std_avg[2,0:6], 'xb-')
axes.plot(array_num_std_avg[2,6], array_mass_std_avg[2,6], 'ob', label = '100K mfa')
axes.grid()
axes.set_xlabel(r"average number coeff. var. $\overline{{\rm CV}(N(t))}$")
axes.set_ylabel(r"average mass coeff. var. $\overline{{\rm CV}(M(t))}$")

axes.annotate(r"$N_{\rm p} = 10^3$", (array_num_std_avg[0,0], array_mass_std_avg[0,0]),
              verticalalignment="bottom", horizontalalignment="center",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$N_{\rm p} = 10^4$", (array_num_std_avg[1,0], array_mass_std_avg[1,0]),
              verticalalignment="bottom", horizontalalignment="center",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$N_{\rm p} = 10^5$", (array_num_std_avg[2,0], array_mass_std_avg[2,0]),
              verticalalignment="bottom", horizontalalignment="center",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')

axes.annotate(r"$\alpha = +1$", (array_num_std_avg[0,0], array_mass_std_avg[0,0]),
              verticalalignment="center", horizontalalignment="left",
              xytext=(5, 0), textcoords='offset points',
              bbox = dict(facecolor='white', edgecolor='white'))
axes.annotate(r"$\alpha = 0$", (array_num_std_avg[0,1], array_mass_std_avg[0,1]),
              verticalalignment="center", horizontalalignment="left",
              xytext=(5, 0), textcoords='offset points',
              bbox = dict(facecolor='white', edgecolor='white'))
axes.annotate(r"$\alpha = -1$", (array_num_std_avg[0,2], array_mass_std_avg[0,2]),
              verticalalignment="center", horizontalalignment="left",
              xytext=(5, 5), textcoords='offset points',
              bbox = dict(facecolor='white', edgecolor='white'))
axes.annotate(r"$\alpha = -2$", (array_num_std_avg[0,3], array_mass_std_avg[0,3]),
              verticalalignment="top", horizontalalignment="center",
              xytext=(0, -5), textcoords='offset points',
              bbox = dict(facecolor='white', edgecolor='white'))
axes.annotate(r"$\alpha = -3$", (array_num_std_avg[0,4], array_mass_std_avg[0,4]),
              verticalalignment="bottom", horizontalalignment="center",
              xytext=(-4, 20), textcoords='offset points',
              bbox = dict(facecolor='white', edgecolor='white'),
              arrowprops=dict(arrowstyle='->'))
axes.annotate(r"$\alpha = -4$", (array_num_std_avg[0,5], array_mass_std_avg[0,5]),
              verticalalignment="bottom", horizontalalignment="center",
              xytext=(-5, 25), textcoords='offset points',
              bbox = dict(facecolor='white', edgecolor='white'),
              arrowprops=dict(arrowstyle='->'))

axes.annotate(r"MFA", (array_num_std_avg[0,6], array_mass_std_avg[0,6]),
              verticalalignment="center", horizontalalignment="left",
              xytext=(5, -2), textcoords='offset points',
              bbox = dict(facecolor='white', edgecolor='white'))

figure.savefig("figs/cv_mass_cv_num.pdf")

