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

    diam = numpy.loadtxt(os.path.join(dirname, "diam.txt"))
    sect_num = numpy.loadtxt(os.path.join(dirname, "num_sect.txt"))
    sect_1_num = numpy.loadtxt(os.path.join(dirname, "num_1_sect.txt"))
    sect_2_num = numpy.loadtxt(os.path.join(dirname, "num_2_sect.txt"))
    part_num = numpy.loadtxt(os.path.join(dirname, "num_dist_mean.txt"))
    part_1_num = numpy.loadtxt(os.path.join(dirname, "num_1_dist_mean.txt"))
    part_1_num_ci = numpy.loadtxt(os.path.join(dirname, "num_1_dist_ci.txt"))
    part_2_num = numpy.loadtxt(os.path.join(dirname, "num_2_dist_mean.txt"))
    part_2_num_ci = numpy.loadtxt(os.path.join(dirname, "num_2_dist_ci.txt"))

    axes.plot(diam * 1e6, sect_1_num, 'b-')
    axes.errorbar(diam * 1e6, part_1_num, error_bars(part_1_num, part_1_num_ci), fmt='b.')
    axes.plot(diam * 1e6, sect_2_num, 'r-')
    axes.errorbar(diam * 1e6, part_2_num, error_bars(part_2_num, part_2_num_ci), fmt='r.')

    mpl_helper.label_plot_line(axes, diam * 1e6, sect_1_num,
                               0.23, r"population 1",
                               horizontalalignment="left",
                               verticalalignment="bottom",
                               draw_box=False)
    mpl_helper.label_plot_line(axes, diam * 1e6, sect_2_num,
                               0.047, r"population 2",
                               horizontalalignment="left",
                               verticalalignment="top",
                               draw_box=False)

    axes.set_xscale('log')
    axes.set_yscale('log')
    axes.set_xlim(left=1e-2, right=1e0)
    axes.set_ylim(bottom=1e4, top=1e12)
    if do_xlabel:
        axes.set_xlabel(r'diameter $D\ /\ {\rm \mu m}$')
    axes.set_ylabel(r'$dN/d\ln D\ /\ {\rm m^{-3}}$')
    axes.grid(True)
    #axes.set_yticks([1e-4, 1e-2, 1e0, 1e2, 1e4, 1e6, 1e8, 1e10])
    mpl_helper.axes_boxed_text(axes, title, "upper left")


if not os.path.exists(config.fig_dirname):
    os.mkdir(config.fig_dirname)

(figure, axes_array) = mpl_helper.make_fig_array(2, 1, figure_width=4, vert_sep=0.2)

run_name = "1k_flat"
axes = axes_array[1][0]
dirname = os.path.join(config.run_dirname, run_name)
make_plot(axes, dirname, "equal weight", do_xlabel=False)

run_name = "1k_flat_source"
axes = axes_array[0][0]
dirname = os.path.join(config.run_dirname, run_name)
make_plot(axes, dirname, "equal number", do_xlabel=True)

mpl_helper.remove_fig_array_axes(axes_array)
figure.savefig(os.path.join(config.fig_dirname, "multi_size_dist_comb.pdf"))
