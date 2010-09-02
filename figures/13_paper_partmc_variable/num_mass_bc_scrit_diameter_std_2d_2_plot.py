#!/usr/bin/env python

import os, sys
import scipy.io
import numpy as np

sys.path.append("../../tool")
import partmc
import mpl_helper
import matplotlib
import config

x_array_bc = np.loadtxt("data/2d_bc_10K_wei-1_12_x_values.txt") * 1e6
y_array_bc = np.loadtxt("data/2d_bc_10K_wei-1_12_y_values.txt") * 100

x_array_scrit = np.loadtxt("data/2d_scrit_10K_wei-1_12_x_values.txt") * 1e6
y_array_scrit = np.loadtxt("data/2d_scrit_10K_wei-1_12_y_values.txt")

num_bc_std = np.loadtxt("data/2d_bc_10K_wei-1_12_hist_std_norm_num.txt")
num_scrit_std = np.loadtxt("data/2d_scrit_10K_wei-1_12_hist_std_norm_num.txt")

num_bc_std = np.ma.masked_invalid(num_bc_std)
num_scrit_std = np.ma.masked_invalid(num_scrit_std)

(figure, axes_array, cbar_axes_array) \
    = mpl_helper.make_fig_array(1,2, figure_width=config.figure_width_double, 
                                left_margin=0.5,right_margin=0.2,top_margin=0.8,
                                vert_sep=0.2, horiz_sep=0.7,
                                colorbar="individual", colorbar_location="top", share_y_axes=False)

axes = axes_array[0][0]
cbar_axes = cbar_axes_array[0][0]
p = axes.pcolor(x_array_bc, y_array_bc, num_bc_std.transpose(), 
                norm = matplotlib.colors.LogNorm(vmin=5e-2, vmax=10),linewidths = 0.1)
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylabel(r"$w_{\rm BC} / \%$")
axes.set_ylim(0, 80)
axes.set_yticks([0, 20, 40, 60, 80])
axes.set_xlim(5e-3, 5)
axes.set_xlabel(r"diameter $D$ / $\rm \mu m$")
axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),orientation='horizontal')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"${\rm CV}(n_{\rm BC}(D, w))$ / $???$")

axes = axes_array[0][1]
cbar_axes = cbar_axes_array[0][1]
p = axes.pcolor(x_array_scrit, y_array_scrit, num_scrit_std.transpose(),
                norm = matplotlib.colors.LogNorm(vmin=5e-2, vmax=10), linewidths = 0.1)
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_ylabel(r"$S_{\rm c} / \%$")
axes.set_ylim(1e-3,1e2)
axes.set_xlim(5e-3, 5)
axes.set_xlabel(r"diameter $D$ / $\rm \mu m$")
axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),orientation='horizontal')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"${\rm CV}(n_{\rm  S}(D, S_{\rm c}))$ / $???$")

figure.savefig("figs/num_mass_bc_scrit_diameter_std_2d.pdf")

