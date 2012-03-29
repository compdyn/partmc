#!/usr/bin/env python

import numpy
import sys, os
sys.path.append("../../tool")
import mpl_helper
import matplotlib
import config

colors = ['b', 'r', 'g']
shapes = ['x', '.', '+']

(figure, axes) = mpl_helper.make_fig(figure_width=5, right_margin=1.6, top_margin=0.4)

handles = []
labels = []
saved_shapes = []
shape_handles = []
shape_labels = []
pop_2_improvements = []
pop_1_degradations = []
for (i_part, (n_part, n_part_name, n_part_tex)) in enumerate(config.n_part_list):
    line_x = []
    err_x = []
    line_y = []
    err_y = []
    for (i_ratio, (ratio_type, ratio, ratio_tex)) in enumerate(config.ratio_list):
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

        plot_extra = False
        if ratio == "0.999":
            plot_extra = True
            shape = 'o'
            shape_labels.append("equal weight")
            num_1_err_equal_weight = num_1_err_mean
            num_2_err_equal_weight = num_2_err_mean
        if ratio == "0.5":
            plot_extra = True
            shape = 's'
            shape_labels.append("equal number")
            num_1_err_equal_number = num_1_err_mean
            num_2_err_equal_number = num_2_err_mean
            pop_2_improvement = num_2_err_equal_weight / num_2_err_equal_number
            pop_1_degradation = num_1_err_equal_number / num_1_err_equal_weight
            print "pop 2 improvement factor: %f" % pop_2_improvement
            print "pop 1 degradation factor: %f" % pop_1_degradation
            pop_2_improvements.append(pop_2_improvement)
            pop_1_degradations.append(pop_1_degradation)
        if plot_extra:
            plotline = axes.plot([num_1_err_mean], [num_2_err_mean],
                                 color=colors[i_part], marker=shape, markerfacecolor='w',
                                 markeredgecolor=colors[i_part], markersize=8, markeredgewidth=1,
                                 linestyle='None')
            if shape not in saved_shapes:
                saved_shapes.append(shape)
                shape_handles.append(matplotlib.lines.Line2D([], [],
                                                             color=colors[i_part], marker=shape, markerfacecolor='w',
                                                             markeredgecolor='k', markersize=8, markeredgewidth=1,
                                                             linestyle='None'))
                #shape_handles.append(plotline)

    #(plotline, caplines, barlinecols) = axes.errorbar(line_x, line_y, err_y, err_x,
    #                                                  fmt=(colors[i_part] + ".-"))
    plotline = axes.plot(line_x, line_y, colors[i_part] + ".-")
    handles.append(plotline)
    labels.append("%s" % n_part_tex)

print "mean pop 2 improvement factor: %f" % numpy.array(pop_2_improvements).mean()
print "mean pop 1 degradation factor: %f" % numpy.array(pop_1_degradations).mean()

axes.set_xscale('log')
axes.set_yscale('log')
axes.set_xlabel(r'population 1 error $E[\|n_1 - n_{1, \rm s}\|_2]$')
axes.set_ylabel(r'population 2 error $E[\|n_2 - n_{2, \rm s}\|_2]$')
axes.grid(True)

(ax_x0, ax_y0) = axes.transAxes.transform_point((0, 0))
(ax_x1, ax_y1) = axes.transAxes.transform_point((1, 1))
upper_left_legend = figure.transFigure.inverted().transform_point((ax_x1 + 10, ax_y1))
lower_left_legend = figure.transFigure.inverted().transform_point((ax_x1 + 10, ax_y0))
figure.legend(handles, labels, loc='upper left', numpoints=3, bbox_to_anchor=upper_left_legend, borderaxespad=0)
figure.legend(shape_handles, shape_labels, loc='lower left', numpoints=1, bbox_to_anchor=lower_left_legend, borderaxespad=0)

figure.savefig(os.path.join(config.fig_dirname, "multi_errors.pdf"))
