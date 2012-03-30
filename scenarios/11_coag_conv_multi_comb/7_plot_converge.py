#!/usr/bin/env python

import numpy
import sys, os
sys.path.append("../../tool")
import mpl_helper
import matplotlib
import config

colors = ['b', 'r', 'g', 'm', 'c', 'y', 'k']

(figure, axes) = mpl_helper.make_fig(figure_width=5, right_margin=1.8, top_margin=0.4)

handles = []
labels = []
i_color = 0
for (i_ratio, (ratio_type, ratio, ratio_tex)) in enumerate(reversed(config.ratio_list)):
    if ratio_type not in config.plot_ratio_list:
        continue
    line1_x = []
    line1_y = []
    err1_y = []
    line2_x = []
    line2_y = []
    err2_y = []
    for (i_part, (n_part, n_part_name, n_part_tex)) in enumerate(config.n_part_list):
        weight_type = "nummass_source"
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

        line1_x.append(int(n_part))
        line1_y.append(num_1_err_mean)
        err1_y.append(num_1_err_ci)

        line2_x.append(int(n_part))
        line2_y.append(num_2_err_mean)
        err2_y.append(num_2_err_ci)

    (plotline, caplines, barlinecols) = axes.errorbar(line1_x, line1_y, err1_y,
                                                      fmt=(colors[i_color] + ".-"))
    handles.append(plotline)
    labels.append(r"$a = 1$, " + ratio_tex)

    (plotline, caplines, barlinecols) = axes.errorbar(line2_x, line2_y, err2_y,
                                                      fmt=(colors[i_color] + ".--"))
    handles.append(plotline)
    labels.append(r"$a = 2$, " + ratio_tex)

    i_color += 1

axes.set_xscale('log')
axes.set_yscale('log')
axes.set_xlabel(r'number of particles $N_{\rm p}$')
axes.set_ylabel(r'error $E[\|n - n_{\rm s}\|_2]$')
axes.grid(True)
figure.legend(handles, labels, loc='center right', numpoints=2)

figure.savefig(os.path.join(config.fig_dirname, "multi_converge.pdf"))
