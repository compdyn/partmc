#!/usr/bin/env python

import numpy
import sys, os
sys.path.append("../../tool")
import mpl_helper
import config

single = {}
multi = {}

colors = ['b', 'r', 'g']

for run in config.all_runs():
    dirname = os.path.join(config.run_dirname, run["name"])
    print dirname

    stats_filename = os.path.join(dirname, "stats.txt")
    stats = numpy.loadtxt(stats_filename)

    num_err_mean = stats[0]
    num_err_ci = stats[1]
    mass_err_mean = stats[2]
    mass_err_ci = stats[3]
    num_var = stats[4]
    mass_var = stats[5]

    if run["weight_type"] == "power":
        if run["n_part"] not in single.keys():
            single[run["n_part"]] = ([], [], [], [], [], [])
        single[run["n_part"]][0].append(num_err_mean)
        single[run["n_part"]][1].append(num_err_ci)
        single[run["n_part"]][2].append(mass_err_mean)
        single[run["n_part"]][3].append(mass_err_ci)
        single[run["n_part"]][4].append(num_var)
        single[run["n_part"]][5].append(mass_var)
    elif run["weight_type"] == "nummass":
        multi[run["n_part"]] = (num_err_mean, num_err_ci, mass_err_mean, mass_err_ci, num_var, mass_var)

(figure, axes) = mpl_helper.make_fig(right_margin=1.8)

handles = []
labels = []
for (i, n_part) in enumerate(single.keys()):
    #handles.append(axes.plot(single[n_part][4], single[n_part][5], colors[i] + 'x-'))
    handles.append(axes.plot(single[n_part][0], single[n_part][2], colors[i] + 'x-'))
    axes.errorbar(single[n_part][0], single[n_part][2], fmt=None, ecolor='k',
                  xerr=single[n_part][1], yerr=single[n_part][3])
    labels.append(r'$N = 10^%d$ single' % int(numpy.log10(int(n_part))))

    #handles.append(axes.plot(multi[n_part][4], multi[n_part][5], colors[i] + 'o'))
    handles.append(axes.plot(multi[n_part][0], multi[n_part][2], colors[i] + 'o'))
    axes.errorbar(multi[n_part][0], multi[n_part][2], fmt=None, ecolor='k',
                  xerr=multi[n_part][1], yerr=multi[n_part][3])
    labels.append(r'$N = 10^%d$ multi' % int(numpy.log10(int(n_part))))

axes.set_xscale('log')
axes.set_yscale('log')
axes.set_xlabel(r'mean number error $E[\|n - n_{\rm s}\|_2]$')
axes.set_ylabel(r'mean mass error $E[\|m - m_{\rm s}\|_2]$')
figure.legend(handles, labels, loc='center right')
axes.grid(True)

filename = "brownian_boomerang.pdf"
figure.savefig(filename)
