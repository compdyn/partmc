#!/usr/bin/env python

import os, sys
import numpy as np

sys.path.append("../../tool")
import partmc
import mpl_helper
import matplotlib
import config

top_label_width = 30 # 32.8
top_label_height = 6 # 6.8
top_label_voffset = 13.5 # 15.5

x_array = np.loadtxt("data/1d_brown_12_x_values.txt") * 1e6
num_avg1 = np.loadtxt("data/1d_brown_12_hist_array_gav_num.txt") / 1e6
num_std1 = np.loadtxt("data/1d_brown_12_e_bars_num.txt") / 1e6

mass_avg1 = np.loadtxt("data/1d_brown_12_hist_array_gav_mass.txt") * 1e9
mass_std1 = np.loadtxt("data/1d_brown_12_e_bars_mass.txt") * 1e9

sect_array_num = np.loadtxt("/Users/nriemer/subversion/partmc/trunk/local_scenarios/brownian_test_paper/out/brownian_sect_size_num.txt")
sect_array_mass = np.loadtxt("/Users/nriemer/subversion/partmc/trunk/local_scenarios/brownian_test_paper/out/brownian_sect_size_mass.txt")

(figure, axes_array) = mpl_helper.make_fig_array(1,2, figure_width=config.figure_width_double, 
                                                 top_margin=0.41, bottom_margin=0.45,
                                                 left_margin=1.07, right_margin=0.65,
                                                 horiz_sep=0.3, vert_sep=0.3,
                                                 share_y_axes=False)

axes = axes_array[0][0]
axes.plot(sect_array_num[:,0]*1e6, sect_array_num[:,12]*np.log(10)/1e6, 'b-')
axes.errorbar(x_array, num_avg1, num_std1, fmt=None, ecolor='r')
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylabel(r"number conc. $n(D)$ / $\rm cm^{-3}$")
axes.set_xlim(1e-2, 1)
axes.set_ylim(0, 5e4)
axes.grid(True)
#axes.text(-0.3, 0.5, r'$\alpha = -1$', horizontalalignment='center',
#           verticalalignment='center', transform=axes.transAxes,
#           rotation='vertical', bbox=dict(edgecolor='black', facecolor='white',
#                                          boxstyle="round,pad=0.5"))
axes.text(0.5, 1.15, r'number', horizontalalignment='center',
           verticalalignment='baseline', transform=axes.transAxes, )
#           bbox=dict(edgecolor='red', facecolor='white', alpha=0.8, linewidth=0.1,
#                     boxstyle="round,pad=0.5"))
pos = axes.transAxes.transform_point((0.5, 1)) * 72.0 / figure.get_dpi()
axes.add_artist(matplotlib.patches.FancyBboxPatch((pos[0] - top_label_width / 2.0,
                                                   pos[1] + top_label_voffset),
                                                  top_label_width, top_label_height,
                                                  boxstyle=matplotlib.patches.BoxStyle.Round(pad=3),
                                                  facecolor='white', edgecolor='black',
                                                  transform=None, clip_on=False))

axes = axes_array[0][1]
axes.plot(sect_array_mass[:,0]*1e6, sect_array_mass[:,12]*np.log(10)*1e9, 'b-')
axes.errorbar(x_array, mass_avg1,mass_std1, fmt=None, ecolor='r')
axes.set_xscale("log")
axes.set_yscale("linear")
axes.set_ylabel(r"mass conc. $m(D)$ / $(\rm \upmu g \ m^{-3})$")
axes.set_xlim(1e-2, 1)
axes.set_ylim(0, 1.2e2)
axes.grid(True)
axes.yaxis.tick_right()
axes.yaxis.set_label_position('right')
axes.text(0.5, 1.15, r'mass', horizontalalignment='center',
           verticalalignment='baseline', transform=axes.transAxes)
pos = axes.transAxes.transform_point((0.5, 1)) * 72.0 / figure.get_dpi()
axes.add_artist(matplotlib.patches.FancyBboxPatch((pos[0] - top_label_width / 2.0,
                                                   pos[1] + top_label_voffset),
                                                  top_label_width, top_label_height,
                                                  boxstyle=matplotlib.patches.BoxStyle.Round(pad=3),
                                                  facecolor='white', edgecolor='black',
                                                  transform=None, clip_on=False))

mpl_helper.remove_fig_array_axes(axes_array, remove_y_axes=False)

figure.savefig("figs/num_mass_brown_diameter_1d.pdf")

