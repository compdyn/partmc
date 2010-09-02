#!/usr/bin/env python

import os, sys
import config
import scipy.io
import numpy as np

sys.path.append("../../tool")
import partmc
import mpl_helper
import matplotlib

bc_std_overall = np.zeros([21,4]) # 12 for weighting schemes, 4 for ss-values
ccn_std_overall = np.zeros([21,4])

i_counter = 0
for counter in ["ss1", "ss2", "ss3", "ss4"]:
    f1 = "data/ccn_std_overall_%s.txt" % counter
    f2 = "data/bc_std_overall_%s.txt" % counter
    ccn_std_overall[:,i_counter] = np.loadtxt(f1)
    bc_std_overall[:,i_counter] = np.loadtxt(f2)
    i_counter += 1

(figure, axes_array) = mpl_helper.make_fig_array(2,1, figure_width=config.figure_width_single, vert_sep=0.2,axis_ratio=1)
axes = axes_array[1][0]
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_ylim(1e-3, 10)
axes.set_ylabel(r"$\overline{{\rm CV}(yyy)}$")
axes.grid()
axes.plot(bc_std_overall[0:6,2], ccn_std_overall[0:6,2], 'r-x', label = '1K S_c = 0.1%')
axes.plot(bc_std_overall[6:12,2], ccn_std_overall[6:12,2], 'g-x', label = '10K S_c = 0.1%')
axes.plot(bc_std_overall[12:18,2], ccn_std_overall[12:18,2], 'b-x', label = '100K S_c = 0.1%')

axes = axes_array[0][0]
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_xlabel(r"$\overline{{\rm CV}(xxx)}$")
axes.set_ylabel(r"$\overline{{\rm CV}(yyy)}$")
axes.grid()

axes.plot(bc_std_overall[12:18,1], ccn_std_overall[12:18,1], 'm-x', label = '100K S_c = 0.01%')
axes.plot(bc_std_overall[12:18,2], ccn_std_overall[12:18,2], 'b-x', label = '100K S_c = 0.1%')
axes.plot(bc_std_overall[12:18,3], ccn_std_overall[12:18,3], 'k-x', label = '100K S_c = 0.5%')

#axes.plot(bc_std_overall[0:6,1], ccn_std_overall[0:6,1],'g-x', label = '1K S_c = 0.01%')

#axes.plot(bc_std_overall[0:6,3], ccn_std_overall[0:6,3], 'm-x', label = '1K S_c = 0.5%')

#axes.plot(bc_std_overall[18,1], ccn_std_overall[18,1], 'go', label = '1K S_c = 0.01% mfa')
#axes.plot(bc_std_overall[18,2], ccn_std_overall[18,2], 'bo', label = '1K S_c = 0.1% mfa')
#axes.plot(bc_std_overall[18,3], ccn_std_overall[18,3], 'mo', label = '1K S_c = 0.5% mfa')

#axes.plot(bc_std_overall[6:12,1], ccn_std_overall[6:12,1], 'g--x', label = '10K S_c = 0.01%')

#axes.plot(bc_std_overall[6:12,3], ccn_std_overall[6:12,3], 'm--x', label = '10K S_c = 0.5%')

#axes.plot(bc_std_overall[19,1], ccn_std_overall[19,1], 'go', label = '10K S_c = 0.01% mfa')
#axes.plot(bc_std_overall[19,2], ccn_std_overall[19,2], 'bo', label = '10K S_c = 0.1% mfa')
#axes.plot(bc_std_overall[19,3], ccn_std_overall[19,3], 'mo', label = '10K S_c = 0.5% mfa')


#axes.plot(bc_std_overall[20,1], ccn_std_overall[20,1], 'go', label = '100K S_c = 0.01% mfa')
#axes.plot(bc_std_overall[20,2], ccn_std_overall[20,2], 'bo', label = '100K S_c = 0.1% mfa')
#axes.plot(bc_std_overall[20,3], ccn_std_overall[20,3], 'mo', label = '100K S_c = 0.5% mfa')

#mpl_helper.label_plot_line(axes, bc_std_overall[0,0:6], ccn_std_overall[0,0:6] , 0.03, "1K",
#                                  verticalalignment="bottom", horizontalalignment="right")
#mpl_helper.label_plot_line(axes, bc_std_overall[1,0:6], ccn_std_overall[1,0:6] , 0.01, "10K",
#                                  verticalalignment="bottom", horizontalalignment="right")
#mpl_helper.label_plot_line(axes, bc_std_overall[2,0:6], ccn_std_overall[2,0:6] , 0.002, "100K",
#                                  verticalalignment="bottom", horizontalalignment="right")

#mpl_helper.label_plot_line(axes, bc_std_overall[0,0:6], ccn_std_overall[0,0:6] , 0.03, "1",
#                                  verticalalignment="bottom", horizontalalignment="left")
#mpl_helper.label_plot_line(axes, bc_std_overall[0,0:6], ccn_std_overall[0,0:6] , 0.02, "0",
#                                  verticalalignment="bottom", horizontalalignment="left")
#mpl_helper.label_plot_line(axes, bc_std_overall[0,0:6], ccn_std_overall[0,0:6] , 0.04, "-1",
#                                  verticalalignment="bottom", horizontalalignment="left")
#mpl_helper.label_plot_line(axes, bc_std_overall[0,0:6], ccn_std_overall[0,0:6] , 0.05, "-2",
#                                  verticalalignment="bottom", horizontalalignment="left")
#mpl_helper.label_plot_line(axes, bc_std_overall[0,0:6], ccn_std_overall[0,0:6] , 0.06, "-3",
#                                  verticalalignment="bottom", horizontalalignment="left")
#mpl_helper.label_plot_line(axes, bc_std_overall[0,0:6], ccn_std_overall[0,0:6] , 0.07, "-4",
#                                  verticalalignment="bottom", horizontalalignment="left")

mpl_helper.remove_fig_array_axes(axes_array)
figure.savefig("figs/cv_scav_bc_cv_ccn_0.1.pdf")

