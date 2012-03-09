#!/usr/bin/env python

import numpy
import sys, os
sys.path.append("../../tool")
import mpl_helper
import config

def make_plot(sect_num_file, sect_mass_file, part_num_file, part_mass_file, color, out_file, title):
    (figure, axes_array) \
        = mpl_helper.make_fig_array(2, 1, figure_width=4, vert_sep=0.2, top_margin=0.5)

    ########################################

    sect_num = numpy.loadtxt(sect_num_file)
    part_num = numpy.loadtxt(part_num_file)
    axes = axes_array[1][0]
    axes.plot(sect_num[:,0] * 1e3, sect_num[:,1], 'k--')
    axes.plot(part_num[:,0] * 1e3, part_num[:,1], color + 'x')
    axes.plot(sect_num[:,0] * 1e3, sect_num[:,2], 'k-')
    axes.plot(part_num[:,0] * 1e3, part_num[:,2], color + 'o')

    mpl_helper.label_plot_line(axes, sect_num[:,0] * 1e3, sect_num[:,1],
                               0.1, r"$t = 0\rm\ min$",
                               horizontalalignment="left",
                               verticalalignment="bottom")
    mpl_helper.label_plot_line(axes, sect_num[:,0] * 1e3, sect_num[:,2],
                               0.1, r"$t = 10\rm\ min$",
                               horizontalalignment="right",
                               verticalalignment="top")

    axes.set_xscale('log')
    axes.set_yscale('log')
    axes.set_ylim(bottom=1e-4, top=1e10)
    axes.set_ylabel(r'$dN/d\ln D\ /\ {\rm m^{-3}}$')
    axes.grid(True)
    axes.set_yticks([1e-4, 1e-2, 1e0, 1e2, 1e4, 1e6, 1e8, 1e10])
    #mpl_helper.axes_boxed_text(axes, "num", "upper left")
    #mpl_helper.axes_boxed_text(axes, title, "upper center", offset_y=-25)

    ########################################

    sect_mass = numpy.loadtxt(sect_mass_file)
    part_mass = numpy.loadtxt(part_mass_file)
    axes = axes_array[0][0]
    axes.plot(sect_mass[:,0] * 1e3, sect_mass[:,1], 'k--')
    axes.plot(part_mass[:,0] * 1e3, part_mass[:,1], color + 'x')
    axes.plot(sect_mass[:,0] * 1e3, sect_mass[:,2], 'k-')
    axes.plot(part_mass[:,0] * 1e3, part_mass[:,2], color + 'o')

    mpl_helper.label_plot_line(axes, sect_mass[:,0] * 1e3, sect_mass[:,1],
                               0.1, r"$t = 0\rm\ min$",
                               horizontalalignment="left",
                               verticalalignment="bottom")
    mpl_helper.label_plot_line(axes, sect_mass[:,0] * 1e3, sect_mass[:,2],
                               0.1, r"$t = 10\rm\ min$",
                               horizontalalignment="right",
                               verticalalignment="top")

    axes.set_xscale('log')
    axes.set_yscale('log')
    axes.set_xlim(left=5e-4, right=2e1)
    axes.set_ylim(bottom=1e-12, top=1e0)
    axes.set_xlabel(r'diameter $D\ /\ {\rm mm}$')
    axes.set_ylabel(r'$dM/d\ln D\ /\ ({\rm kg\ m^{-3}})$')
    axes.grid(True)
    axes.set_yticks([1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1e0])
    #mpl_helper.axes_boxed_text(axes, "mass", "upper left")

    ########################################

    mpl_helper.remove_fig_array_axes(axes_array)

    figure.savefig(out_file)

if not os.path.exists(config.fig_dirname):
    os.mkdir(config.fig_dirname)

for run in config.all_runs():
    print run["name"]
    dirname = os.path.join(config.run_dirname, run["name"])
    figname = os.path.join(config.fig_dirname, run["name"] + ".pdf")

    make_plot(os.path.join(dirname, "sect_aero_size_num.txt"),
              os.path.join(dirname, "sect_aero_size_mass.txt"),
              os.path.join(dirname, "part_0001_0001_aero_size_num.txt"),
              os.path.join(dirname, "part_0001_0001_aero_size_mass.txt"),
              "b", figname, run["n_part_tex"] + " particles")


#make_plot("runs/1k_power0/sect_aero_size_num.txt", "runs/1k_power0/sect_aero_size_mass.txt",
#          "runs/1k_power0/part_0001_0001_aero_size_num.txt", "runs/1k_power0/part_0001_0001_aero_size_mass.txt",
#          "g", os.path.join(config.fig_dirname, "sedi_1k_power0.pdf"), r"$10^3$ particles")
#make_plot("out_old_1e6/sedi_sect_size_num.txt", "out_old_1e6/sedi_sect_size_mass.txt",
#          "out_old_1e6/sedi_part_size_num.txt", "out_old_1e6/sedi_part_size_mass.txt",
#          "b", "figs/sedi_old_1e6.pdf", r"$10^6$ particles")
#make_plot("out_new/sedi_sect_aero_size_num.txt", "out_new/sedi_sect_aero_size_mass.txt",
#          "out_new/sedi_part_0001_size_num.txt", "out_new/sedi_part_0001_size_mass.txt",
#          "r", "figs/sedi_new.pdf", r"$10^3$ superparticles")
