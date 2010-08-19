#!/usr/bin/env python2.5

import os, sys
import config
import numpy as np

sys.path.append("../../tool")
import partmc
import mpl_helper
import matplotlib

top_label_width = 32.8
top_label_height = 6.8
top_label_voffset = 15.5

x_array = np.loadtxt("data/1d_10K_wei+1_12_x_values.txt") * 1e6
num_avg1 = np.loadtxt("data/1d_10K_wei+1_12_hist_array_gav_num.txt") / 1e6
num_avg2 = np.loadtxt("data/1d_10K_wei-1_12_hist_array_gav_num.txt") / 1e6
num_avg3 = np.loadtxt("data/1d_10K_wei-4_12_hist_array_gav_num.txt") / 1e6

num_std1 = np.loadtxt("data/1d_10K_wei+1_12_e_bars_num.txt") / 1e6
num_std2 = np.loadtxt("data/1d_10K_wei-1_12_e_bars_num.txt") / 1e6
num_std3 = np.loadtxt("data/1d_10K_wei-4_12_e_bars_num.txt") / 1e6

mass_avg1 = np.loadtxt("data/1d_10K_wei+1_12_hist_array_gav_mass.txt") * 1e9
mass_avg2 = np.loadtxt("data/1d_10K_wei-1_12_hist_array_gav_mass.txt") * 1e9
mass_avg3 = np.loadtxt("data/1d_10K_wei-4_12_hist_array_gav_mass.txt") * 1e9

mass_std1 = np.loadtxt("data/1d_10K_wei+1_12_e_bars_mass.txt") * 1e9
mass_std2 = np.loadtxt("data/1d_10K_wei-1_12_e_bars_mass.txt") * 1e9
mass_std3 = np.loadtxt("data/1d_10K_wei-4_12_e_bars_mass.txt") * 1e9

(figure, axes_array) = mpl_helper.make_fig_array(3,2, figure_width=config.figure_width_double, 
                                                 top_margin=0.41, bottom_margin=0.45,
                                                 left_margin=1.07, right_margin=0.65,
                                                 horiz_sep=0.3, vert_sep=0.3,
                                                 share_y_axes=False)

axes = axes_array[2][0]
axes.errorbar(x_array, num_avg1, num_std1, fmt='b-')
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_ylabel(r"$n(D)$ / $\rm cm^{-3}$")
axes.set_ylim(1e-3, 1e5)
axes.grid(True)
axes.text(-0.4, 0.5, r'$\alpha = +1$', horizontalalignment='center',
           verticalalignment='center', transform=axes.transAxes,
           rotation='vertical', bbox=dict(edgecolor='black', facecolor='white',
                                          boxstyle="round,pad=0.5"))
axes.text(0.5, 1.15, r'number', horizontalalignment='center',
           verticalalignment='baseline', transform=axes.transAxes)
#           bbox=dict(edgecolor='red', facecolor='white', alpha=0.8, linewidth=0.1,
#                     boxstyle="round,pad=0.5"))
pos = axes.transAxes.transform_point((0.5, 1)) * 72.0 / figure.get_dpi()
axes.add_artist(matplotlib.patches.FancyBboxPatch((pos[0] - top_label_width / 2.0,
                                                   pos[1] + top_label_voffset),
                                                  top_label_width, top_label_height,
                                                  boxstyle=matplotlib.patches.BoxStyle.Round(pad=5),
                                                  facecolor='white', edgecolor='black',
                                                  transform=None, clip_on=False))

axes = axes_array[1][0]
axes.errorbar(x_array, num_avg2, num_std2, fmt='b-')
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_ylabel(r"$n(D)$ / $\rm cm^{-3}$")
axes.set_ylim(1e-3, 1e5)
axes.set_xlim(5e-3, 5)
axes.grid(True)
axes.text(-0.4, 0.5, r'$\alpha = -1$', horizontalalignment='center',
           verticalalignment='center', transform=axes.transAxes,
           rotation='vertical', bbox=dict(edgecolor='black', facecolor='white',
                                          boxstyle="round,pad=0.5"))

axes = axes_array[0][0]
axes.errorbar(x_array, num_avg3, num_std3, fmt='b-')
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_ylabel(r"$n(D)$ / $\rm cm^{-3}$")
axes.set_ylim(1e-3, 1e5)
axes.set_xlim(5e-3, 5)
axes.set_xlim(5e-3, 5)
axes.grid(True)
axes.set_xlabel(r"diameter / $\rm \mu m$")
axes.text(-0.4, 0.5, r'$\alpha = -4$', horizontalalignment='center',
           verticalalignment='center', transform=axes.transAxes,
           rotation='vertical', bbox=dict(edgecolor='black', facecolor='white',
                                          boxstyle="round,pad=0.5"))

axes = axes_array[2][1]
axes.errorbar(x_array, mass_avg1,mass_std1, fmt='r-')
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_ylabel(r"$m(D)$ / $\rm \mu g \, m^{-3}$")
axes.set_ylim(1e-4, 1e2)
axes.grid(True)
axes.yaxis.tick_right()
axes.yaxis.set_label_position('right')
axes.text(0.5, 1.15, r'mass', horizontalalignment='center',
           verticalalignment='baseline', transform=axes.transAxes)
pos = axes.transAxes.transform_point((0.5, 1)) * 72.0 / figure.get_dpi()
axes.add_artist(matplotlib.patches.FancyBboxPatch((pos[0] - top_label_width / 2.0,
                                                   pos[1] + top_label_voffset),
                                                  top_label_width, top_label_height,
                                                  boxstyle=matplotlib.patches.BoxStyle.Round(pad=5),
                                                  facecolor='white', edgecolor='black',
                                                  transform=None, clip_on=False))
axes = axes_array[1][1]
axes.errorbar(x_array, mass_avg2, mass_std2, fmt='r-')
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_ylabel(r"$m(D)$ / $\rm \mu g \, m^{-3}$")
axes.set_ylim(1e-4, 1e2)
axes.grid(True)
axes.yaxis.tick_right()
axes.yaxis.set_label_position('right')

axes = axes_array[0][1]
axes.errorbar(x_array, mass_avg3, mass_std3, fmt='r-')
axes.set_xscale("log")
axes.set_yscale("log")
axes.set_ylabel(r"$m(D)$ / $\rm \mu g \, m^{-3}$")
axes.set_ylim(1e-4, 1e2)
axes.set_xlim(5e-3, 5)
axes.grid(True)
axes.yaxis.tick_right()
axes.yaxis.set_label_position('right')
axes.set_xlabel(r"diameter / $\rm \mu m$")

mpl_helper.remove_fig_array_axes(axes_array, remove_y_axes=False)

figure.savefig("figs/num_mass_diameter_1d.pdf")

