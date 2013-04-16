#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import scipy.io
import numpy as np
import matplotlib.pyplot
from matplotlib.patches import Polygon
from matplotlib.patches import FancyArrowPatch

matplotlib.rc('text.latex', preamble='\usepackage{rotating}')

(figure, axes) = mpl_helper.make_fig(figure_width=4, axis_ratio=1)

ncf = scipy.io.netcdf_file("urban_plume/out/urban_plume2_process.nc")
avg_part_entropy1 = np.exp(ncf.variables["avg_part_entropy"].data)
entropy_of_avg_part1 = np.exp(ncf.variables["entropy_of_avg_part"].data)
ncf.close()

axes.plot((avg_part_entropy1-1)/9, (entropy_of_avg_part1-1)/9, "b-", linewidth=2)

axes.set_xlabel(r"$(D_{\alpha}-1)/(A-1)$")
axes.set_xlim([0.2,1.])
axes.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0])
axes.set_ylim([0.2,1])
axes.set_yticks([0.2, 0.4, 0.6, 0.8, 1])
axes.set_ylabel(r"$(D_{\gamma}-1)/(A-1)$")

axes.annotate("",
            xy=((avg_part_entropy1[30]-1)/9,(entropy_of_avg_part1[30]-1)/9), xycoords='data',
            xytext=((avg_part_entropy1[25]-1)/9,(entropy_of_avg_part1[25]-1)/9), textcoords='data',
            arrowprops=dict(arrowstyle="simple, head_width=0.5, head_length=1, tail_width=0.001",
                            connectionstyle="arc3", ec ="b", fc = "b", ),
            )


x_values = np.linspace(0,2,21)

axes.plot(x_values, x_values/0.5, "k", linewidth=0.5)
axes.plot(x_values, x_values/0.75, "k", linewidth=0.5)
axes.plot(x_values, x_values, "k", linewidth=0.5)

axes.annotate(r"$\chi=0.5$", (0.45, 0.9),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')

axes.annotate(r"$\chi=0.75$", (0.68, 0.9),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')

axes.annotate(r"$\chi=1$", (0.9,0.9),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')

axes.grid(True)
figure.savefig("d_g_versus_d_a_urban_plume_poster.pdf")


