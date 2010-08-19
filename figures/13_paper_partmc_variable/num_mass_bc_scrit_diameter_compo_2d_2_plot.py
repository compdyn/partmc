#!/usr/bin/env python

import os, sys
import config
import scipy.io
import numpy as np

sys.path.append("../../tool")
import partmc
import mpl_helper
import matplotlib

x_array_bc = np.loadtxt("data/2d_bc_compo_12_x_values.txt") * 1e6
y_array_bc = np.loadtxt("data/2d_bc_compo_12_y_values.txt") * 100

x_array_scrit = np.loadtxt("data/2d_scrit_compo_12_x_values.txt") * 1e6
y_array_scrit = np.loadtxt("data/2d_scrit_compo_12_y_values.txt") 

num_bc = np.loadtxt("data/2d_bc_compo_12_average_num.txt") / 1e6
num_scrit = np.loadtxt("data/2d_scrit_compo_12_average_num.txt") / 1e6

mass_bc = np.loadtxt("data/2d_bc_compo_12_average_mass.txt") * 1e9
mass_scrit = np.loadtxt("data/2d_scrit_compo_12_average_mass.txt") * 1e9

(figure, axes_array, cbar_axes_array) = mpl_helper.make_fig_array(2,2, figure_width=config.figure_width_double, 
                                                                  left_margin=0.7, right_margin=1.6,
                                                                  vert_sep=0.3, horiz_sep=0.8,
                                                                  share_y_axes=False, colorbar=True)

axes = axes_array[1][0]
cbar_axes = cbar_axes_array[1][0]
p = axes.pcolor(x_array_bc, y_array_bc, num_bc.transpose(),linewidths = 0.1, norm=matplotlib.colors.LogNorm())
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylabel(r"BC mass frac. $w_{\rm BC}$ / \%")
axes.set_ylim(0, 80)
axes.grid(True)
figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext())
cbar_axes.set_ylabel(r"$n(D,w)$ / $\rm cm^{-3}$")

axes = axes_array[0][0]
cbar_axes = cbar_axes_array[0][0]
p = axes.pcolor(x_array_bc, y_array_bc, mass_bc.transpose(),linewidths = 0.1, norm=matplotlib.colors.LogNorm())
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylabel(r"BC mass frac. $w_{\rm BC}$ / \%")
axes.set_ylim(0, 80)
axes.set_xlim(5e-3, 5)
axes.set_xlabel(r"diameter $D$ / $\rm \mu m$")
axes.grid(True)
figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext())
cbar_axes.set_ylabel(r"$m(D,w)$ / $\rm \mu g \ m^{-3}$")

axes = axes_array[1][1]
cbar_axes = cbar_axes_array[1][1]
p = axes.pcolor(x_array_scrit, y_array_scrit, num_scrit.transpose(),linewidths = 0.1, norm=matplotlib.colors.LogNorm())
axes.set_xscale("log")
axes.set_yscale("log")
axes.yaxis.tick_right()
axes.yaxis.set_label_position('right')
axes.set_ylabel(r"crit. supersat. $S_{\rm crit}$ / \%")
axes.set_ylim(1e-3, 1e2)
axes.grid(True)
figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext())
cbar_axes.set_ylabel(r"$n(D,S)$ / $\rm cm^{-3}$")

axes = axes_array[0][1]
cbar_axes = cbar_axes_array[0][1]
p = axes.pcolor(x_array_scrit, y_array_scrit, mass_scrit.transpose(),linewidths = 0.1, norm=matplotlib.colors.LogNorm())
axes.set_xscale("log")
axes.set_yscale("log")
axes.yaxis.tick_right()
axes.yaxis.set_label_position('right')
axes.set_ylabel(r"crit. supersat. $S_{\rm crit}$ / \%")
axes.set_ylim(1e-3, 1e2)
axes.set_xlim(5e-3, 5)
axes.set_xlabel(r"diameter $D$ / $\rm \mu m$")
axes.grid(True)
figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext())
cbar_axes.set_ylabel(r"$m(D,S)$ / $\rm \mu g \ m^{-3}$")

mpl_helper.remove_fig_array_axes(axes_array, remove_y_axes=False)

figure.savefig("figs/num_mass_bc_scrit_diameter_compo_2d.pdf")

