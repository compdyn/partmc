#!/usr/bin/env python

import numpy
import sys, os
sys.path.append("../../tool")
import mpl_helper
import config

def make_plot(dirname, out_file, title):
    (figure, axes_array) \
        = mpl_helper.make_fig_array(2, 1, figure_width=4, vert_sep=0.2, top_margin=0.5)

    ########################################

    diam = numpy.loadtxt(os.path.join(dirname, "diam.txt"))
    sect_num = numpy.loadtxt(os.path.join(dirname, "num_sect.txt"))
    sect_1_num = numpy.loadtxt(os.path.join(dirname, "num_1_sect.txt"))
    sect_2_num = numpy.loadtxt(os.path.join(dirname, "num_2_sect.txt"))
    part_num = numpy.loadtxt(os.path.join(dirname, "num_dist_mean.txt"))
    part_1_num = numpy.loadtxt(os.path.join(dirname, "num_1_dist_mean.txt"))
    part_2_num = numpy.loadtxt(os.path.join(dirname, "num_2_dist_mean.txt"))
    axes = axes_array[1][0]
    axes.plot(diam * 1e6, sect_num, 'k--')
    axes.plot(diam * 1e6, part_num, 'k+')
    axes.plot(diam * 1e6, sect_1_num + sect_2_num, 'g:')
    axes.plot(diam * 1e6, sect_1_num, 'b:')
    axes.plot(diam * 1e6, part_1_num, 'bx')
    axes.plot(diam * 1e6, sect_2_num, 'r--')
    axes.plot(diam * 1e6, part_2_num, 'rx')

    #mpl_helper.label_plot_line(axes, sect_num[:,0] * 1e3, sect_num[:,1],
    #                           0.1, r"$t = 0\rm\ min$",
    #                           horizontalalignment="left",
    #                           verticalalignment="bottom")
    #mpl_helper.label_plot_line(axes, sect_num[:,0] * 1e3, sect_num[:,2],
    #                           0.1, r"$t = 10\rm\ min$",
    #                           horizontalalignment="right",
    #                           verticalalignment="top")

    axes.set_xscale('log')
    axes.set_yscale('log')
    axes.set_ylim(bottom=1e8, top=1e12)
    axes.set_ylabel(r'$dN/d\ln D\ /\ {\rm m^{-3}}$')
    axes.grid(True)
    #axes.set_yticks([1e-4, 1e-2, 1e0, 1e2, 1e4, 1e6, 1e8, 1e10])
    mpl_helper.axes_boxed_text(axes, "num", "upper left")
    mpl_helper.axes_boxed_text(axes, title, "upper center", offset_y=-25)

    ########################################

    diam = numpy.loadtxt(os.path.join(dirname, "diam.txt"))
    sect_mass = numpy.loadtxt(os.path.join(dirname, "mass_sect.txt"))
    sect_1_mass = numpy.loadtxt(os.path.join(dirname, "mass_1_sect.txt"))
    sect_2_mass = numpy.loadtxt(os.path.join(dirname, "mass_2_sect.txt"))
    part_mass = numpy.loadtxt(os.path.join(dirname, "mass_dist_mean.txt"))
    part_1_mass = numpy.loadtxt(os.path.join(dirname, "mass_1_dist_mean.txt"))
    part_2_mass = numpy.loadtxt(os.path.join(dirname, "mass_2_dist_mean.txt"))
    axes = axes_array[0][0]
    axes.plot(diam * 1e6, sect_mass, 'k--')
    axes.plot(diam * 1e6, part_mass, 'k+')
    axes.plot(diam * 1e6, sect_1_mass + sect_2_mass, 'g:')
    axes.plot(diam * 1e6, sect_1_mass, 'b:')
    axes.plot(diam * 1e6, part_1_mass, 'bx')
    axes.plot(diam * 1e6, sect_2_mass, 'r--')
    axes.plot(diam * 1e6, part_2_mass, 'rx')

    #mpl_helper.label_plot_line(axes, sect_mass[:,0] * 1e3, sect_mass[:,1],
    #                           0.1, r"$t = 0\rm\ min$",
    #                           horizontalalignment="left",
    #                           verticalalignment="bottom")
    #mpl_helper.label_plot_line(axes, sect_mass[:,0] * 1e3, sect_mass[:,2],
    #                           0.1, r"$t = 10\rm\ min$",
    #                           horizontalalignment="right",
    #                           verticalalignment="top")

    axes.set_xscale('log')
    axes.set_yscale('log')
    axes.set_xlim(left=1e-2, right=1e0)
    axes.set_ylim(bottom=1e-12, top=1e-6)
    axes.set_xlabel(r'diameter $D\ /\ {\rm \mu m}$')
    axes.set_ylabel(r'$dM/d\ln D\ /\ ({\rm kg\ m^{-3}})$')
    axes.grid(True)
    #axes.set_yticks([1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e0])
    mpl_helper.axes_boxed_text(axes, "mass", "upper left")

    ########################################

    mpl_helper.remove_fig_array_axes(axes_array)

    figure.savefig(out_file)

if not os.path.exists(config.fig_dirname):
    os.mkdir(config.fig_dirname)

for run in config.all_runs():
    print run["name"]
    dirname = os.path.join(config.run_dirname, run["name"])
    figname = os.path.join(config.fig_dirname, run["name"] + ".pdf")

    make_plot(dirname, figname, run["n_part_tex"] + " particles")
