#!/usr/bin/env python

import numpy
import sys, os
sys.path.append("../../tool")
import mpl_helper
import config

def error_bars(mean, ci):
    scale = 1 + ci / mean
    return numpy.vstack((mean - mean / scale, ci))

def make_plot(axes1, axes2, dirname, title1, title2, do_xlabel, top_label):

    diam = numpy.loadtxt(os.path.join(dirname, "diam.txt"))
    sect_num = numpy.loadtxt(os.path.join(dirname, "num_sect.txt"))
    sect_1_num = numpy.loadtxt(os.path.join(dirname, "num_1_sect.txt"))
    sect_2_num = numpy.loadtxt(os.path.join(dirname, "num_2_sect.txt"))
    part_num = numpy.loadtxt(os.path.join(dirname, "num_dist_mean.txt"))
    part_1_num = numpy.loadtxt(os.path.join(dirname, "num_1_dist_mean.txt"))
    part_1_num_ci = numpy.loadtxt(os.path.join(dirname, "num_1_dist_ci.txt"))
    part_2_num = numpy.loadtxt(os.path.join(dirname, "num_2_dist_mean.txt"))
    part_2_num_ci = numpy.loadtxt(os.path.join(dirname, "num_2_dist_ci.txt"))

    sect_mass = numpy.loadtxt(os.path.join(dirname, "mass_sect.txt"))
    sect_1_mass = numpy.loadtxt(os.path.join(dirname, "mass_1_sect.txt"))
    sect_2_mass = numpy.loadtxt(os.path.join(dirname, "mass_2_sect.txt"))
    part_mass = numpy.loadtxt(os.path.join(dirname, "mass_dist_mean.txt"))
    part_1_mass = numpy.loadtxt(os.path.join(dirname, "mass_1_dist_mean.txt"))
    part_1_mass_ci = numpy.loadtxt(os.path.join(dirname, "mass_1_dist_ci.txt"))
    part_2_mass = numpy.loadtxt(os.path.join(dirname, "mass_2_dist_mean.txt"))
    part_2_mass_ci = numpy.loadtxt(os.path.join(dirname, "mass_2_dist_ci.txt"))

    axes1.plot(diam * 1e6, sect_1_num, 'b-')
    axes1.errorbar(diam * 1e6, part_1_num, error_bars(part_1_num, part_1_num_ci), fmt='b.')
    axes1.plot(diam * 1e6, sect_2_num, 'r-')
    axes1.errorbar(diam * 1e6, part_2_num, error_bars(part_2_num, part_2_num_ci), fmt='r.')

    axes2.plot(diam * 1e6, sect_1_mass, 'b-')
    axes2.errorbar(diam * 1e6, part_1_mass, error_bars(part_1_mass, part_1_mass_ci), fmt='b.')
    axes2.plot(diam * 1e6, sect_2_mass, 'r-')
    axes2.errorbar(diam * 1e6, part_2_mass, error_bars(part_2_mass, part_2_mass_ci), fmt='r.')

    #mpl_helper.label_plot_line(axes1, diam * 1e6, sect_1_num,
    #                           0.23, r"population 1",
    #                           horizontalalignment="left",
    #                           verticalalignment="bottom",
    #                           draw_box=False)
    #mpl_helper.label_plot_line(axes1, diam * 1e6, sect_2_num,
    #                           0.047, r"population 2",
    #                           horizontalalignment="left",
    #                           verticalalignment="top",
    #                           draw_box=False)
    #mpl_helper.label_plot_line(axes2, diam * 1e6, sect_1_mass,
    #                           0.23, r"population 1",
    #                           horizontalalignment="left",
    #                           verticalalignment="bottom",
    #                           draw_box=False)
    #mpl_helper.label_plot_line(axes2, diam * 1e6, sect_2_mass,
    #                           0.047, r"population 2",
    #                           horizontalalignment="left",
    #                           verticalalignment="top",
    #                           draw_box=False)

    axes1.set_xscale('log')
    axes1.set_yscale('log')
    axes2.set_xscale('log')
    axes2.set_yscale('log')
    axes1.set_xlim(left=1e-2, right=1e1)
    axes2.set_xlim(left=1e-2, right=1e1)
    axes1.set_ylim(bottom=1e4, top=1e12)
    axes2.set_ylim(bottom=1e-16, top=1e-6)
    if do_xlabel:
        axes1.set_xlabel(r'diameter $D\ /\ {\rm \mu m}$')
        axes2.set_xlabel(r'diameter $D\ /\ {\rm \mu m}$')
    axes1.set_ylabel(r'$n$ / ${\rm m^{-3}}$')
    axes2.yaxis.tick_right()
    axes2.yaxis.set_label_position('right')
    axes2.set_ylabel(r'$m$ / $({\rm kg\ m^{-3}})$')
    axes1.grid(True)
    axes2.grid(True)
    axes1.set_yticks([1e4, 1e6, 1e8, 1e10, 1e12])
    axes2.set_yticks([1e-16, 1e-14, 1e-12, 1e-10, 1e-8, 1e-6])
    mpl_helper.axes_boxed_text(axes1, title1, "upper left", box_padding=0.3, offset_x=7, offset_y=7)
    mpl_helper.axes_boxed_text(axes2, title1, "upper left", box_padding=0.3, offset_x=7, offset_y=7)
    mpl_helper.axes_boxed_text(axes1, title2, "upper right", box_padding=0.3, offset_x=7, offset_y=7)
    mpl_helper.axes_boxed_text(axes2, title2, "upper right", box_padding=0.3, offset_x=7, offset_y=7)
    if top_label:
        mpl_helper.axes_boxed_text(axes1, "number", "upper center", offset_y=-20)
        mpl_helper.axes_boxed_text(axes2, "mass", "upper center", offset_y=-20)

if not os.path.exists(config.fig_dirname):
    os.mkdir(config.fig_dirname)

(figure, axes_array) = mpl_helper.make_fig_array(4, 2, figure_width=5, share_y_axes=False,
                                                 vert_sep=0.2, horiz_sep=0.2,
                                                 right_margin=0.7, top_margin=0.4, left_margin=0.65)

run_name = "1k_flat"
axes1 = axes_array[3][0]
axes2 = axes_array[3][1]
dirname = os.path.join(config.run_dirname, run_name)
make_plot(axes1, axes2, dirname, r"$\alpha = 0$", r"$r = 1$", do_xlabel=False, top_label=True)

run_name = "1k_w1-w2_flat_source"
axes1 = axes_array[2][0]
axes2 = axes_array[2][1]
dirname = os.path.join(config.run_dirname, run_name)
make_plot(axes1, axes2, dirname, r"$\alpha = 0$", r"$r = 10^3$", do_xlabel=False, top_label=False)

run_name = "1k_w1-1e-3_nummass"
axes1 = axes_array[1][0]
axes2 = axes_array[1][1]
dirname = os.path.join(config.run_dirname, run_name)
make_plot(axes1, axes2, dirname, r"comb.", r"$r = 1$", do_xlabel=False, top_label=False)

run_name = "1k_w1-w2_nummass_source"
axes1 = axes_array[0][0]
axes2 = axes_array[0][1]
dirname = os.path.join(config.run_dirname, run_name)
make_plot(axes1, axes2, dirname, r"comb.", r"$r = 10^3$", do_xlabel=True, top_label=False)

mpl_helper.remove_fig_array_axes(axes_array, remove_y_axes=False)
figure.savefig(os.path.join(config.fig_dirname, "both_size_dist.pdf"))
