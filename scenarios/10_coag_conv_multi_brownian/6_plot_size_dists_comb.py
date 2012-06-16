#!/usr/bin/env python

import numpy
import sys, os
sys.path.append("../../tool")
import mpl_helper
import config

def error_bars(mean, ci):
    scale = 1 + ci / mean
    return numpy.vstack((mean - mean / scale, ci))

def make_plot(axes, dirname, title, do_xlabel):

    diam = numpy.loadtxt(os.path.join(dirname, "diam.txt")) * 1e6
    sect_num = numpy.loadtxt(os.path.join(dirname, "num_sect.txt")) / 1e6
    sect_1_num = numpy.loadtxt(os.path.join(dirname, "num_1_sect.txt")) / 1e6
    sect_2_num = numpy.loadtxt(os.path.join(dirname, "num_2_sect.txt")) / 1e6
    part_num = numpy.loadtxt(os.path.join(dirname, "num_dist_mean.txt")) / 1e6
    part_1_num = numpy.loadtxt(os.path.join(dirname, "num_1_dist_mean.txt")) / 1e6
    part_1_num_ci = numpy.loadtxt(os.path.join(dirname, "num_1_dist_ci.txt")) / 1e6
    part_2_num = numpy.loadtxt(os.path.join(dirname, "num_2_dist_mean.txt")) / 1e6
    part_2_num_ci = numpy.loadtxt(os.path.join(dirname, "num_2_dist_ci.txt")) / 1e6

    axes.plot(diam, sect_1_num, 'b-')
    axes.plot(diam, sect_2_num, 'r-')

    mpl_helper.label_plot_line(axes, diam, sect_1_num,
                               0.14, r"sub-population 1",
                               horizontalalignment="left",
                               verticalalignment="bottom",
                               draw_box=True,
                               default_offset_mag=8)
    mpl_helper.label_plot_line(axes, diam, sect_2_num,
                               0.045, r"sub-population 2",
                               horizontalalignment="left",
                               verticalalignment="top",
                               draw_box=True,
                               default_offset_mag=8)

    axes.errorbar(diam, part_1_num, error_bars(part_1_num, part_1_num_ci), fmt='b.')
    axes.errorbar(diam, part_2_num, error_bars(part_2_num, part_2_num_ci), fmt='r.')

    axes.set_xscale('log')
    axes.set_yscale('log')
    axes.set_xlim(left=1e-2, right=1e0)
    axes.set_ylim(bottom=1e-2, top=1e6)
    if do_xlabel:
        axes.set_xlabel(r'diameter $D\ /\ {\rm \mu m}$')
    axes.set_ylabel(r'$dN/d\ln D\ /\ {\rm cm^{-3}}$')
    axes.grid(True)
    #axes.set_yticks([1e-4, 1e-2, 1e0, 1e2, 1e4, 1e6, 1e8, 1e10])
    mpl_helper.axes_boxed_text(axes, title, "upper left", box_padding=0.3, offset_x=7, offset_y=7)


if not os.path.exists(config.fig_dirname):
    os.mkdir(config.fig_dirname)

(figure, axes_array) = mpl_helper.make_fig_array(2, 1, figure_width=4, vert_sep=0.2)

run_name = "1k_flat"
axes = axes_array[1][0]
dirname = os.path.join(config.run_dirname, run_name)
make_plot(axes, dirname, r"equal weight, $r = 1$", do_xlabel=False)

run_name = "1k_flat_source"
axes = axes_array[0][0]
dirname = os.path.join(config.run_dirname, run_name)
make_plot(axes, dirname, r"equal number, $r = 10^3$", do_xlabel=True)

mpl_helper.remove_fig_array_axes(axes_array)
figure.savefig(os.path.join(config.fig_dirname, "multi_size_dist_comb.pdf"))
