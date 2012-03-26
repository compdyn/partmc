#!/usr/bin/env python

import numpy
import sys, os
sys.path.append("../../tool")
import mpl_helper
import config

fig_name = os.path.join(config.fig_dirname, "comb_size_dists.pdf")

if not os.path.exists(config.fig_dirname):
    os.mkdir(config.fig_dirname)

n_rows = len(config.plot_weight_list)
(figure, axes_array) \
    = mpl_helper.make_fig_array(n_rows, 2, figure_width=5, share_y_axes=False,
                                vert_sep=0.2, horiz_sep=0.2,
                                right_margin=0.7, top_margin=0.4, left_margin=0.65)

for (i_weight, (weight_type, exponent)) in enumerate(config.plot_weight_list):
    for (i_part, (n_part, n_part_name, n_part_tex)) in enumerate(config.n_part_list):
        if n_part != "1000":
            continue
        if weight_type == "power":
            name = "%s_%s%s" % (n_part_name, weight_type, exponent)
        else:
            name = "%s_%s" % (n_part_name, weight_type)

        print name

        i_row = n_rows - i_weight - 1
        if i_row == 0:
            alpha_text = "comb."
        else:
            alpha_text = r"$\alpha = %s$" % exponent

        dirname = os.path.join(config.run_dirname, name)
        sect_num_file = os.path.join(dirname, "sect_aero_size_num.txt")
        sect_mass_file = os.path.join(dirname, "sect_aero_size_mass.txt")
        part_num_file = os.path.join(dirname, "part_0001_0001_aero_size_num.txt")
        part_mass_file = os.path.join(dirname, "part_0001_0001_aero_size_mass.txt")
    
        ########################################
    
        sect_num = numpy.loadtxt(sect_num_file)
        part_num = numpy.loadtxt(part_num_file)
        axes = axes_array[i_row][0]
        axes.plot(sect_num[:,0] * 1e3, sect_num[:,1], 'b-')
        axes.plot(part_num[:,0] * 1e3, part_num[:,1], 'b.')
        axes.plot(sect_num[:,0] * 1e3, sect_num[:,2], 'r-')
        axes.plot(part_num[:,0] * 1e3, part_num[:,2], 'r.')
    
        mpl_helper.label_plot_line(axes, sect_num[:,0] * 1e3, sect_num[:,1],
                                   0.15, r"$t = 0$",
                                   horizontalalignment="left",
                                   verticalalignment="bottom")
        mpl_helper.label_plot_line(axes, sect_num[:,0] * 1e3, sect_num[:,2],
                                   0.1, r"$t = 10$",
                                   horizontalalignment="right",
                                   verticalalignment="top")
    
        axes.set_xscale('log')
        axes.set_yscale('log')
        axes.set_xlim(left=5e-4, right=2e1)
        axes.set_ylim(bottom=1e-4, top=1e10)
        if i_row == 0:
            axes.set_xlabel(r'diameter $D\ /\ {\rm mm}$')
        axes.set_ylabel(r'$n$ / ${\rm m^{-3}}$')
        axes.grid(True)
        axes.set_yticks([1e-4, 1e-2, 1e0, 1e2, 1e4, 1e6, 1e8, 1e10])
        mpl_helper.axes_boxed_text(axes, alpha_text, "upper right",
                                   box_padding=0.2, offset_x=5, offset_y=5)
        if i_weight == 0:
            mpl_helper.axes_boxed_text(axes, "number", "upper center", offset_y=-20)
        
        ########################################
        
        sect_mass = numpy.loadtxt(sect_mass_file)
        part_mass = numpy.loadtxt(part_mass_file)
        axes = axes_array[i_row][1]
        axes.plot(sect_mass[:,0] * 1e3, sect_mass[:,1], 'b-')
        axes.plot(part_mass[:,0] * 1e3, part_mass[:,1], 'b.')
        axes.plot(sect_mass[:,0] * 1e3, sect_mass[:,2], 'r-')
        axes.plot(part_mass[:,0] * 1e3, part_mass[:,2], 'r.')
    
        mpl_helper.label_plot_line(axes, sect_mass[:,0] * 1e3, sect_mass[:,1],
                                   0.1, r"$t = 0$",
                                   horizontalalignment="left",
                                   verticalalignment="bottom")
        mpl_helper.label_plot_line(axes, sect_mass[:,0] * 1e3, sect_mass[:,2],
                                   0.1, r"$t = 10$",
                                   horizontalalignment="right",
                                   verticalalignment="top")
    
        axes.set_xscale('log')
        axes.set_yscale('log')
        axes.set_xlim(left=5e-4, right=2e1)
        axes.set_ylim(bottom=1e-12, top=1e0)
        if i_row == 0:
            axes.set_xlabel(r'diameter $D\ /\ {\rm mm}$')
        axes.yaxis.tick_right()
        axes.yaxis.set_label_position('right')
        axes.set_ylabel(r'$m$ / $({\rm kg\ m^{-3}})$')
        axes.grid(True)
        axes.set_yticks([1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e0])
        mpl_helper.axes_boxed_text(axes, alpha_text, "upper left",
                                   box_padding=0.2, offset_x=5, offset_y=5)
        if i_weight == 0:
            mpl_helper.axes_boxed_text(axes, "mass", "upper center", offset_y=-20)

mpl_helper.remove_fig_array_axes(axes_array, remove_y_axes=False)
figure.savefig(fig_name)
