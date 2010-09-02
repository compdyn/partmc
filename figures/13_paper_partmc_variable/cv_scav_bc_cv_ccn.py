#!/usr/bin/env python

import os, sys
import config
import scipy.io
import numpy as np

sys.path.append("../../tool")
import partmc
import mpl_helper
import matplotlib

bc_std_overall = np.zeros([21,4]) # 21 for weighting schemes, 4 for ss-values
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
axes.set_ylabel(r"$\overline{{\rm CV}(M_{\rm BC}(S_{\rm env}))}$")
axes.grid()
axes.plot(ccn_std_overall[0:6,2], bc_std_overall[0:6,2], 'r-x', label = '1K S_c = 0.1%')
axes.plot(ccn_std_overall[6:12,2], bc_std_overall[6:12,2], 'g-x', label = '10K S_c = 0.1%')
axes.plot(ccn_std_overall[12:18,2], bc_std_overall[12:18,2], 'b-x', label = '100K S_c = 0.1%')
axes.plot(ccn_std_overall[18,2], bc_std_overall[18,2], 'ro', label = '1K S_c = 0.1% mfa')
axes.plot(ccn_std_overall[19,2], bc_std_overall[19,2], 'go', label = '10K S_c = 0.1% mfa')
axes.plot(ccn_std_overall[20,2], bc_std_overall[20,2], 'bo', label = '100K S_c = 0.1% mfa')

mpl_helper.label_plot_line(axes, ccn_std_overall[0:6,2], bc_std_overall[0:6,2] , 0.4, r"$N_{\rm p} = 10^3$",
                                  verticalalignment="center", horizontalalignment="left")
mpl_helper.label_plot_line(axes, ccn_std_overall[6:12,2], bc_std_overall[6:12,2] , 0.2, r"$N_{\rm p} = 10^4$",
                                  verticalalignment="center", horizontalalignment="left")
mpl_helper.label_plot_line(axes, ccn_std_overall[12:18,2], bc_std_overall[12:18,2] , 0.1, r"$N_{\rm p} =10^5$",
                                  verticalalignment="center", horizontalalignment="left")
mpl_helper.label_plot_line(axes, ccn_std_overall[12:18,2], bc_std_overall[12:18,2] , 0.02, r"$\alpha = 1$",
                                  verticalalignment="bottom", horizontalalignment="right")
mpl_helper.label_plot_line(axes, ccn_std_overall[12:18,2], bc_std_overall[12:18,2] , 0.009, r"$\alpha = 0$",
                                  verticalalignment="bottom", horizontalalignment="right")
mpl_helper.label_plot_line(axes, ccn_std_overall[12:18,2], bc_std_overall[12:18,2] , 0.006, r"$\alpha = -1$",
                                  verticalalignment="center", horizontalalignment="right")
mpl_helper.label_plot_line(axes, ccn_std_overall[12:18,2], bc_std_overall[12:18,2] , 0.01, r"$\alpha = -2$",
                                  verticalalignment="top", horizontalalignment="right")
mpl_helper.label_plot_line(axes, ccn_std_overall[12:18,2], bc_std_overall[12:18,2] , 0.03, r"$\alpha = -3$",
                                  verticalalignment="top", horizontalalignment="left")
mpl_helper.label_plot_line(axes, ccn_std_overall[12:18,2], bc_std_overall[12:18,2] , 0.01, r"$\alpha = -4$",
                                  verticalalignment="bottom", horizontalalignment="left")

axes = axes_array[0][0]
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_xlabel(r"$\overline{{\rm CV}(N_{\rm CCN}(S_{\rm env}))}$")
axes.set_ylabel(r"$\overline{{\rm CV}(M_{\rm BC}(S_{\rm env}))}$")
axes.grid()

axes.plot(ccn_std_overall[12:18,1], bc_std_overall[12:18,1], 'm-x', label = '100K S_c = 0.01%')
axes.plot(ccn_std_overall[12:18,2], bc_std_overall[12:18,2], 'b-x', label = '100K S_c = 0.1%')
axes.plot(ccn_std_overall[12:18,3], bc_std_overall[12:18,3], 'k-x', label = '100K S_c = 0.5%')

axes.plot(ccn_std_overall[20,1], bc_std_overall[20,1], 'mo', label = '100K S_c = 0.01% mfa')
axes.plot(ccn_std_overall[20,2], bc_std_overall[20,2], 'bo', label = '100K S_c = 0.1% mfa')
axes.plot(ccn_std_overall[20,3], bc_std_overall[20,3], 'ko', label = '100K S_c = 0.5% mfa')

mpl_helper.label_plot_line(axes, ccn_std_overall[12:18,1], bc_std_overall[12:18,1] , 0.3, r"$S_{\rm c} = 0.01$ %",
                                  verticalalignment="center", horizontalalignment="left")
mpl_helper.label_plot_line(axes, ccn_std_overall[12:18,2], bc_std_overall[12:18,2] , 0.2, r"$S_{\rm c} = 0.1$ %",
                                  verticalalignment="center", horizontalalignment="left")
mpl_helper.label_plot_line(axes, ccn_std_overall[12:18,3], bc_std_overall[12:18,3] , 0.1, r"$S_{\rm c} =0.5$ %",
                           verticalalignment="center", horizontalalignment="left")
mpl_helper.remove_fig_array_axes(axes_array)
figure.savefig("figs/cv_scav_bc_cv_ccn.pdf")

