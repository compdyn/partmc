#!/usr/bin/env python

import numpy
import sys, os
sys.path.append("../../tool")
import mpl_helper
import config

colors = ['b', 'r', 'g']
shapes = ['x', 'o', '+']

(figure, axes) = mpl_helper.make_fig(right_margin=2)

for run in config.all_runs():
    dirname = os.path.join(config.run_dirname, run["name"])
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

    axes.plot([num_1_err_mean], [num_2_err_mean], colors[run["i_part"]] + shapes[run["i_weight"]],
              label="%s, %s" % (run["n_part_tex"], run["weight_type"].replace("_", "-")))

axes.set_xscale('log')
axes.set_yscale('log')
axes.set_xlabel(r'mean number 1 error $E[\|n_1 - n_{1, \rm s}\|_2]$')
axes.set_ylabel(r'mean number 2 error $E[\|n_2 - n_{2, \rm s}\|_2]$')
(handles, labels) = axes.get_legend_handles_labels()
figure.legend(handles, labels, loc='center right', numpoints=1)
axes.grid(True)

figure.savefig(os.path.join(config.fig_dirname, "errors.pdf"))
