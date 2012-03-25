#!/usr/bin/env python

import numpy
import sys, os
sys.path.append("../../tool")
import mpl_helper
import config

colors = ['b', 'r', 'g']
shapes = ['x', '.', '+']

(figure, axes) = mpl_helper.make_fig(right_margin=2)

handles = []
labels = []
for (i_part, (n_part, n_part_name, n_part_tex)) in enumerate(config.n_part_list):
    line_x = []
    err_x = []
    line_y = []
    err_y = []
    for (i_ratio, (ratio_type, ratio)) in enumerate(config.ratio_list):
        weight_type = "flat_source"
        name = "%s_%s_%s" % (n_part_name, ratio_type, weight_type)
        dirname = os.path.join(config.run_dirname, name)
        print dirname

        stats_filename = os.path.join(dirname, "stats.txt")
        stats = numpy.loadtxt(stats_filename)

        num_1_err_mean = stats[0]
        num_1_err_ci = stats[1]
        num_2_err_mean = stats[2]
        num_2_err_ci = stats[3]
        mass_1_err_mean = stats[4]
        mass_1_err_ci = stats[5]
        mass_2_err_mean = stats[6]
        mass_2_err_ci = stats[7]

        line_x.append(num_1_err_mean)
        line_y.append(num_2_err_mean)

        err_x.append(num_1_err_ci)
        err_y.append(num_2_err_ci)

    (plotline, caplines, barlinecols) = axes.errorbar(line_x, line_y, err_y, err_x,
                                                      fmt=(colors[i_part]))
    handles.append(plotline)
    labels.append("%s" % n_part_tex)

axes.set_xscale('log')
axes.set_yscale('log')
axes.set_xlabel(r'mean number 1 error $E[\|n_1 - n_{1, \rm s}\|_2]$')
axes.set_ylabel(r'mean number 2 error $E[\|n_2 - n_{2, \rm s}\|_2]$')
figure.legend(handles, labels, loc='center right', numpoints=1)
axes.grid(True)

figure.savefig(os.path.join(config.fig_dirname, "errors.pdf"))
