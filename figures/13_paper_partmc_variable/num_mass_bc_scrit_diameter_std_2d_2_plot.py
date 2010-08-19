#!/usr/bin/env python2.5

import os, sys
import config
import scipy.io
import numpy as np

sys.path.append("../../tool")
import partmc
import mpl_helper
import matplotlib

x_array_bc = np.loadtxt("data/2d_bc_10K_wei+1_12_x_values.txt") * 1e6
y_array_bc = np.loadtxt("data/2d_bc_10K_wei+1_12_y_values.txt") * 100

x_array_scrit = np.loadtxt("data/2d_scrit_10K_wei+1_12_x_values.txt") * 1e6
y_array_scrit = np.loadtxt("data/2d_scrit_10K_wei+1_12_y_values.txt")

def min_pos(a):
    ma = np.ma.masked_less_equal(a, 0)
    return ma.min()

num_bc_std = np.loadtxt("data/2d_bc_10K_wei+1_12_hist_std_norm_num.txt")
num_scrit_std = np.loadtxt("data/2d_scrit_10K_wei+1_12_hist_std_norm_num.txt")
print "range: ", min_pos(num_bc_std), num_bc_std.max()
print "range: ", min_pos(num_scrit_std), num_scrit_std.max()

(figure, axes_array, cbar_axes_array) \
    = mpl_helper.make_fig_array(2,1, figure_width=config.figure_width_single, 
                                colorbar=True,
                                left_margin=0.7, right_margin=0.6, vert_sep=0.3)

axes = axes_array[1][0]
cbar_axes = cbar_axes_array[0]
p = axes.pcolor(x_array_bc, y_array_bc, num_bc_std.transpose(), norm = matplotlib.colors.LogNorm(vmin=1e-2, vmax=10),linewidths = 0.1)
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylabel(r"$w_{\rm BC} / \%$")
axes.set_ylim(0, 80)
axes.set_xlim(5e-3, 5)
axes.grid(True)
#figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
#                orientation='vertical')

axes = axes_array[0][0]
axes.pcolor(x_array_scrit, y_array_scrit, num_scrit_std.transpose(),linewidths = 0.1)
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_ylabel(r"$S_{\rm crit} / \%$")
axes.set_ylim(1e-3,1e2)
axes.set_xlim(5e-3, 5)
axes.grid(True)

mpl_helper.remove_fig_array_axes(axes_array)

figure.savefig("figs/num_mass_bc_scrit_diameter_std_2d.pdf")

