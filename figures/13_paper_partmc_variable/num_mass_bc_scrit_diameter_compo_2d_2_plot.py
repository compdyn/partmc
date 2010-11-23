#!/usr/bin/env python

import os, sys
import scipy.io
import numpy as np

sys.path.append("../../tool")
import partmc
import mpl_helper
import matplotlib
import config

num_sc_min = 1e-3
num_sc_max = 4e5
mass_sc_min = 1e-6
mass_sc_max = 4e1

num_bc_min = 4e0
num_bc_max = 4e5
mass_bc_min = 1e-5
mass_bc_max = 4e1

x_array_bc = np.loadtxt("data/2d_bc_compo_12_x_values.txt") * 1e6
y_array_bc = np.loadtxt("data/2d_bc_compo_12_y_values.txt") * 100

x_array_scrit = np.loadtxt("data/2d_scrit_compo_12_x_values.txt") * 1e6
y_array_scrit = np.loadtxt("data/2d_scrit_compo_12_y_values.txt") 

num_bc = np.loadtxt("data/2d_bc_compo_12_average_num.txt") / 1e6
num_scrit = np.loadtxt("data/2d_scrit_compo_12_average_num.txt") / 1e6

mass_bc = np.loadtxt("data/2d_bc_compo_12_average_mass.txt") * 1e9
mass_scrit = np.loadtxt("data/2d_scrit_compo_12_average_mass.txt") * 1e9

(figure, axes_array, cbar_axes_array) = mpl_helper.make_fig_array(2,2, figure_width=config.figure_width_double, 
                                                                  left_margin=0.7, right_margin=0.1,
                                                                  top_margin=0.8, vert_sep=1.5, horiz_sep=0.3,
                                                                  colorbar="individual", colorbar_location="top")

axes = axes_array[1][0]
cbar_axes = cbar_axes_array[1][0]
p = axes.pcolor(x_array_bc, y_array_bc, num_bc.transpose(),linewidths = 0.1, 
                norm=matplotlib.colors.LogNorm(vmin=num_bc_min, vmax=num_bc_max))
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylabel(r"BC mass frac. $w_{\rm BC}$ / \%")
axes.set_ylim(0, 80)
axes.set_xlabel(r"dry diameter $D_{\rm dry}$ / $\rm \upmu m$")
axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                       orientation='horizontal')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"number conc. $n(D_{\rm dry},w_{\rm BC})$ / $\rm cm^{-3}$")
cbar.set_ticks([1e-3, 1e-1, 1e1, 1e3, 1e5])

axes = axes_array[0][0]
cbar_axes = cbar_axes_array[0][0]
p = axes.pcolor(x_array_scrit, y_array_scrit, num_scrit.transpose(),linewidths = 0.1, 
                norm=matplotlib.colors.LogNorm(vmin=num_sc_min, vmax=num_sc_max))
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_ylabel(r"crit. supersat. $S_{\rm c}$ / \%")
axes.set_ylim(1e-3, 1e2)
axes.set_xlim(5e-3, 5)
axes.set_xlabel(r"dry diameter $D_{\rm dry}$ / $\rm \upmu m$")
axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                       orientation='horizontal')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"number conc. $n(D_{\rm dry},S_{\rm c})$ / $\rm cm^{-3}$")
#cbar.set_ticks([1e-3, 1e-1, 1e1, 1e3, 1e5])

axes = axes_array[1][1]
cbar_axes = cbar_axes_array[1][1]
p = axes.pcolor(x_array_bc, y_array_bc, mass_bc.transpose(),linewidths = 0.1, 
                norm=matplotlib.colors.LogNorm(vmin=mass_bc_min, vmax=mass_bc_max))
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_xlabel(r"dry diameter $D_{\rm dry}$ / $\rm \upmu m$")
axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                       orientation='horizontal')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"mass conc. $m^{\rm BC}(D_{\rm dry},w_{\rm BC})$ / $(\rm \upmu g \ m^{-3})$")
cbar.set_ticks([1e-7, 1e-5, 1e-3, 1e-1, 1e1])

axes = axes_array[0][1]
cbar_axes = cbar_axes_array[0][1]
p = axes.pcolor(x_array_scrit, y_array_scrit, mass_scrit.transpose(),linewidths = 0.1, 
                norm=matplotlib.colors.LogNorm(vmin=mass_sc_min, vmax=mass_sc_max))
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_xlim(5e-3, 5)
axes.set_xlabel(r"dry diameter $D_{\rm dry}$ / $\rm \upmu m$")
axes.grid(True)
cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                       orientation='horizontal')
cbar_axes.xaxis.set_label_position('top')
cbar.set_label(r"mass conc. $m^{\rm BC}(D_{\rm dry},S_{\rm c})$ / $(\rm \upmu g \ m^{-3})$")
cbar.set_ticks([1e-7, 1e-5, 1e-3, 1e-1, 1e1])

mpl_helper.remove_fig_array_axes(axes_array, remove_x_axes=False)

figure.savefig("figs/num_mass_bc_scrit_diameter_compo_2d.pdf")

