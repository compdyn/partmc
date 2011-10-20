#!/usr/bin/env python

import numpy
import sys, os
sys.path.append("../../tool")
import mpl_helper
import config

single = {}
multi = {}

for run in config.all_runs():
    dirname = os.path.join(config.run_dirname, run["name"])
    print dirname

    stats_filename = os.path.join(dirname, "stats.txt")
    stats = numpy.loadtxt(stats_filename)

    num_err_mean = stats[0]
    num_err_rel_ci = stats[1]
    mass_err_mean = stats[2]
    mass_err_rel_ci = stats[3]

    if run["weight_type"] == "power":
        if run["n_part"] not in single.keys():
            single[run["n_part"]] = ([], [], [], [])
        single[run["n_part"]][0].append(num_err_mean)
        single[run["n_part"]][1].append(num_err_rel_ci)
        single[run["n_part"]][2].append(mass_err_mean)
        single[run["n_part"]][3].append(mass_err_rel_ci)
    elif run["weight_type"] == "nummass":
        multi[run["n_part"]] = (num_err_mean, num_err_rel_ci, mass_err_mean, mass_err_rel_ci)

(figure, axes) = mpl_helper.make_fig(right_margin=1.8)

handles = []
labels = []
for n_part in single.keys():
    handles.append(axes.plot(single[n_part][0], single[n_part][2], 'bx-'))
    labels.append(r'$N = %s$ single' % n_part)

    handles.append(axes.plot(multi[n_part][0], multi[n_part][2], 'bo'))
    labels.append(r'$N = %s$ multi' % n_part)

axes.set_xscale('log')
axes.set_yscale('log')
axes.set_xlabel(r'mean number error $E[\|n - n_{\rm s}\|_2]$')
axes.set_ylabel(r'mean mass error $E[\|m - m_{\rm s}\|_2]$')
figure.legend(handles, labels, loc='center right')

filename = "boomerang.pdf"
figure.savefig(filename)
