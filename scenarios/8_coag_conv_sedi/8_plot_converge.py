#!/usr/bin/env python

import numpy
import sys, os
sys.path.append("../../tool")
import mpl_helper
import matplotlib
import config

colors = ['b', 'r', 'g', 'm', 'c', 'y', 'k']

(figure, axes) = mpl_helper.make_fig(figure_width=5, right_margin=2.1, top_margin=0.1, left_margin=0.6)
axes2 = axes.twinx()

handles = []
labels = []
i_color = 0
for (weight_type, exponent) in config.plot_weight_list:
    num_line_x = []
    num_line_y = []
    num_err_y = []
    mass_line_x = []
    mass_line_y = []
    mass_err_y = []
    for (i_part, (n_part, n_part_name, n_part_tex)) in enumerate(config.n_part_list):
        if weight_type == "power":
            name = "%s_%s%s" % (n_part_name, weight_type, exponent)
            exponent_tex = r"$\alpha = %s$" % exponent
        else:
            name = "%s_%s" % (n_part_name, weight_type)
            exponent_tex = r"comb."
        dirname = os.path.join(config.run_dirname, name)
        print dirname

        stats_filename = os.path.join(dirname, "stats.txt")
        stats = numpy.loadtxt(stats_filename)

        num_err_mean = stats[0]
        num_err_ci = stats[1]
        mass_err_mean = stats[2]
        mass_err_ci = stats[3]
        num_dist_var = stats[4]
        mass_dist_var = stats[5]

        num_line_x.append(int(n_part))
        num_line_y.append(num_err_mean)
        num_err_y.append(num_err_ci)

        mass_line_x.append(int(n_part))
        mass_line_y.append(mass_err_mean)
        mass_err_y.append(mass_err_ci)

    (plotline, caplines, barlinecols) = axes.errorbar(num_line_x, num_line_y, num_err_y,
                                                      fmt=(colors[i_color] + ".-"))
    handles.append(plotline)
    labels.append(exponent_tex)

    (plotline, caplines, barlinecols) = axes2.errorbar(mass_line_x, mass_line_y, mass_err_y,
                                                      fmt=(colors[i_color] + ".--"))

    i_color += 1

axes.set_xscale('log')
axes.set_yscale('log')
axes2.set_yscale('log')
axes.set_xlabel(r'number of particles $N_{\rm p}$')
axes.set_ylabel(r'error $E[\|n - n_{\rm fv}\|_2]$')
axes2.set_ylabel(r'error $E[\|m - m_{\rm fv}\|_2]$')
axes.grid(True)

line_handles = [
    matplotlib.lines.Line2D([], [], color='k', marker='None', linestyle='-'),
    matplotlib.lines.Line2D([], [], color='k', marker='None', linestyle='--'),
    ]
line_labels = [
    "number",
    "mass",
    ]

(ax_x0, ax_y0) = axes.transAxes.transform_point((0, 0))
(ax_x1, ax_y1) = axes.transAxes.transform_point((1, 1))
upper_left_legend = figure.transFigure.inverted().transform_point((ax_x1 + 70, ax_y1))
lower_left_legend = figure.transFigure.inverted().transform_point((ax_x1 + 70, ax_y0))
figure.legend(handles, labels, loc='upper left', numpoints=2, bbox_to_anchor=upper_left_legend, borderaxespad=0)
figure.legend(line_handles, line_labels, loc='lower left', numpoints=2, bbox_to_anchor=lower_left_legend, borderaxespad=0,
              handlelength=2.5)

figure.savefig(os.path.join(config.fig_dirname, "comb_converge.pdf"))
