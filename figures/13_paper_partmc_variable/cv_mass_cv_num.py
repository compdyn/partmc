#!/usr/bin/env python2.5

import os, sys
import config
import scipy.io
import numpy as np

sys.path.append("../../tool")
import partmc
import mpl_helper
import matplotlib

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

print 'array_num_std_avg ', array_num_std_avg
print 'mfa ', array_num_std_avg[0,6], array_mass_std_avg[0,6]
x_array = [1, 0, -1, -2, -3, -4]

(figure, axes) = mpl_helper.make_fig(figure_width=config.figure_width_single)

axes.set_xscale("log")
axes.set_yscale("log")
axes.plot(array_num_std_avg[0,0:6], array_mass_std_avg[0,0:6], 'xr-', label = '1K')
#axes.plot(array_num_std_avg[0,6], array_mass_std_avg[0,6], 'or', label = '1K mfa')
axes.plot(array_num_std_avg[1,0:6], array_mass_std_avg[1,0:6], 'xg-', label = '10K')
#axes.plot(array_num_std_avg[1,6], array_mass_std_avg[1,6], 'og', label = '10K mfa')
axes.plot(array_num_std_avg[2,0:6], array_mass_std_avg[2,0:6], 'xb-', label = '100K')
#axes.plot(array_num_std_avg[2,6], array_mass_std_avg[2,6], 'ob', label = '100K mfa')
axes.grid()
axes.set_xlabel(r"${\rm CV}(\overline{N(t)})$")
axes.set_ylabel(r"${\rm CV}(\overline{M(t)})$")


mpl_helper.label_plot_line(axes, array_num_std_avg[0,0:6], array_mass_std_avg[0,0:6] , 0.03, "1K",
                                  verticalalignment="bottom", horizontalalignment="right")
mpl_helper.label_plot_line(axes, array_num_std_avg[1,0:6], array_mass_std_avg[1,0:6] , 0.01, "10K",
                                  verticalalignment="bottom", horizontalalignment="right")
mpl_helper.label_plot_line(axes, array_num_std_avg[2,0:6], array_mass_std_avg[2,0:6] , 0.002, "100K",
                                  verticalalignment="bottom", horizontalalignment="right")

mpl_helper.label_plot_line(axes, array_num_std_avg[0,0:6], array_mass_std_avg[0,0:6] , 0.03, "1",
                                  verticalalignment="bottom", horizontalalignment="left")
mpl_helper.label_plot_line(axes, array_num_std_avg[0,0:6], array_mass_std_avg[0,0:6] , 0.02, "0",
                                  verticalalignment="bottom", horizontalalignment="left")
mpl_helper.label_plot_line(axes, array_num_std_avg[0,0:6], array_mass_std_avg[0,0:6] , 0.04, "-1",
                                  verticalalignment="bottom", horizontalalignment="left")
mpl_helper.label_plot_line(axes, array_num_std_avg[0,0:6], array_mass_std_avg[0,0:6] , 0.05, "-2",
                                  verticalalignment="bottom", horizontalalignment="left")
mpl_helper.label_plot_line(axes, array_num_std_avg[0,0:6], array_mass_std_avg[0,0:6] , 0.06, "-3",
                                  verticalalignment="bottom", horizontalalignment="left")
mpl_helper.label_plot_line(axes, array_num_std_avg[0,0:6], array_mass_std_avg[0,0:6] , 0.07, "-4",
                                  verticalalignment="bottom", horizontalalignment="left")

figure.savefig("figs/cv_mass_cv_num.pdf")

