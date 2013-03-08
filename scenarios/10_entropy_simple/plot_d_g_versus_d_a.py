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

(figure, axes) = mpl_helper.make_fig(figure_width=5, axis_ratio=1)

ncf = scipy.io.netcdf_file("emit/out/urban_plume_process.nc")
avg_part_entropy1 = np.exp(ncf.variables["avg_part_entropy"].data)
entropy_of_avg_part1 = np.exp(ncf.variables["entropy_of_avg_part"].data)
ncf.close()

ncf = scipy.io.netcdf_file("emit_mixed/out/urban_plume_process.nc")
avg_part_entropy2 = np.exp(ncf.variables["avg_part_entropy"].data)
entropy_of_avg_part2 = np.exp(ncf.variables["entropy_of_avg_part"].data)
ncf.close()

ncf = scipy.io.netcdf_file("emit_pure/out/urban_plume_process.nc")
avg_part_entropy3 = np.exp(ncf.variables["avg_part_entropy"].data)
entropy_of_avg_part3 = np.exp(ncf.variables["entropy_of_avg_part"].data)
ncf.close()

ncf = scipy.io.netcdf_file("coag/out/urban_plume_process.nc")
avg_part_entropy4 = np.exp(ncf.variables["avg_part_entropy"].data)
entropy_of_avg_part4 = np.exp(ncf.variables["entropy_of_avg_part"].data)
ncf.close()

ncf = scipy.io.netcdf_file("cond_mono/out/urban_plume_process.nc")
avg_part_entropy5 = np.exp(ncf.variables["avg_part_entropy"].data)
entropy_of_avg_part5 = np.exp(ncf.variables["entropy_of_avg_part"].data)
ncf.close()

ncf = scipy.io.netcdf_file("cond_mono_decrease/out/urban_plume_process.nc")
avg_part_entropy6 = np.exp(ncf.variables["avg_part_entropy"].data)
entropy_of_avg_part6 = np.exp(ncf.variables["entropy_of_avg_part"].data)
ncf.close()

ncf = scipy.io.netcdf_file("cond_bidisp/out/urban_plume_process.nc")
avg_part_entropy7 = np.exp(ncf.variables["avg_part_entropy"].data)
entropy_of_avg_part7 = np.exp(ncf.variables["entropy_of_avg_part"].data)
ncf.close()

ncf = scipy.io.netcdf_file("cond_bidisp_decrease/out/urban_plume_process.nc")
avg_part_entropy8 = np.exp(ncf.variables["avg_part_entropy"].data)
entropy_of_avg_part8 = np.exp(ncf.variables["entropy_of_avg_part"].data)
ncf.close()

axes.plot((avg_part_entropy1-1)/2, (entropy_of_avg_part1-1)/2, "b-", linewidth=2)
axes.plot((avg_part_entropy2-1)/3, (entropy_of_avg_part2-1)/3, "b--", linewidth=2)
axes.plot((avg_part_entropy3-1)/2, (entropy_of_avg_part3-1)/2, "c--", linewidth=2)
axes.plot((avg_part_entropy4-1)/1, (entropy_of_avg_part4-1)/1, "r-", linewidth=2)
axes.plot((avg_part_entropy5-1)/3, (entropy_of_avg_part5-1)/3, "k-", linewidth=2)
axes.plot((avg_part_entropy6-1)/3, (entropy_of_avg_part6-1)/3, "k--", linewidth=2)
axes.plot((avg_part_entropy7-1)/3, (entropy_of_avg_part7-1)/3, "g-", linewidth=2)
axes.plot((avg_part_entropy8-1)/3, (entropy_of_avg_part8-1)/3, "g--", linewidth=2)

x_values = np.linspace(0,2,21)

axes.plot(x_values, x_values*1, "k", linewidth=0.5)
axes.plot(x_values, x_values*2, "k", linewidth=0.5)

#axes.annotate(r"1", (avg_part_entropy1[5], entropy_of_avg_part1[5]),
#              verticalalignment="bottom", horizontalalignment="right",
#              bbox = dict(facecolor='white', edgecolor='white'),
#              xytext=(0, 5), textcoords='offset points')


#axes.annotate(r"$\gamma=0$", (x_values[1], 1),
#              verticalalignment="bottom", horizontalalignment="right",
#              bbox = dict(facecolor='white', edgecolor='white'),
#              xytext=(0, 5), textcoords='offset points')

#axes.annotate(r"$\gamma=0.5$", (x_values[5], 2*x_values[5]),
#              verticalalignment="bottom", horizontalalignment="right",
#              bbox = dict(facecolor='white', edgecolor='white'),
#              xytext=(0, 5), textcoords='offset points')

#axes.annotate(r"$\gamma=1$", (x_values[10], x_values[10]),
#              verticalalignment="bottom", horizontalalignment="right",
#              bbox = dict(facecolor='white', edgecolor='white'),
#              xytext=(0, 5), textcoords='offset points')


axes.set_xlabel(r"$(D_{\alpha}-1)/(A-1)$")
axes.set_xlim([0,1.1])
axes.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
axes.set_ylim([0,1.1])
axes.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
axes.set_ylabel(r"$(D_{\gamma}-1)/(A-1)$")
pts = np.array([[0,0], [1,1], [0,1]])
p = Polygon(pts,closed=True,alpha=0.05)
axes.add_patch(p)

axes.annotate("",
            xy=((avg_part_entropy1[45]-1)/2,(entropy_of_avg_part1[45]-1)/2), xycoords='data',
            xytext=((avg_part_entropy1[25]-1)/2,(entropy_of_avg_part1[25]-1)/2), textcoords='data',
            arrowprops=dict(arrowstyle="simple, head_width=0.5, head_length=1, tail_width=0.001",
                            connectionstyle="arc3", ec ="b", fc = "b", ),
            )

axes.annotate("",
            xy=((avg_part_entropy2[6]-1)/3,(entropy_of_avg_part2[6]-1)/3), xycoords='data',
            xytext=((avg_part_entropy2[4]-1)/3,(entropy_of_avg_part2[4]-1)/3), textcoords='data',
            arrowprops=dict(arrowstyle="simple, head_width=0.5, head_length=1, tail_width=0.001",
                            connectionstyle="arc3", ec ="b", fc = "b", ),
            )

axes.annotate("",
            xy=((avg_part_entropy3[10]-1)/2,(entropy_of_avg_part3[10]-1)/2), xycoords='data',
            xytext=((avg_part_entropy3[5]-1)/2,(entropy_of_avg_part3[5]-1)/2), textcoords='data',
            arrowprops=dict(arrowstyle="simple, head_width=0.5, head_length=1, tail_width=0.001",
                            connectionstyle="arc3", ec ="c", fc = "c", ),
            )

axes.annotate("",
            xy=((avg_part_entropy4[15]-1)/1,(entropy_of_avg_part4[15]-1)/1), xycoords='data',
            xytext=((avg_part_entropy4[5]-1)/1,(entropy_of_avg_part4[5]-1)/1), textcoords='data',
            arrowprops=dict(arrowstyle="simple, head_width=0.5, head_length=1, tail_width=0.001",
                            connectionstyle="arc3", ec ="r", fc = "r", ),
            )

axes.annotate("",
            xy=((avg_part_entropy5[7]-1)/3,(entropy_of_avg_part5[7]-1)/3), xycoords='data',
            xytext=((avg_part_entropy5[3]-1)/3,(entropy_of_avg_part5[3]-1)/3), textcoords='data',
            arrowprops=dict(arrowstyle="simple, head_width=0.5, head_length=1, tail_width=0.001",
                            connectionstyle="arc3", ec ="k", fc = "k", ),
            )

axes.annotate("",
            xy=((avg_part_entropy6[8]-1)/3,(entropy_of_avg_part6[8]-1)/3), xycoords='data',
            xytext=((avg_part_entropy6[5]-1)/3,(entropy_of_avg_part6[5]-1)/3), textcoords='data',
            arrowprops=dict(arrowstyle="simple, head_width=0.5, head_length=1, tail_width=0.001",
                            connectionstyle="arc3", ec ="k", fc = "k", ),
            )

axes.annotate("",
            xy=((avg_part_entropy7[16]-1)/3,(entropy_of_avg_part7[16]-1)/3), xycoords='data',
            xytext=((avg_part_entropy7[6]-1)/3,(entropy_of_avg_part7[6]-1)/3), textcoords='data',
            arrowprops=dict(arrowstyle="simple, head_width=0.5, head_length=1, tail_width=0.001",
                            connectionstyle="arc3", ec ="g", fc = "g", ),
            )

axes.annotate("",
            xy=((avg_part_entropy8[15]-1)/3,(entropy_of_avg_part8[15]-1)/3), xycoords='data',
            xytext=((avg_part_entropy8[5]-1)/3,(entropy_of_avg_part8[5]-1)/3), textcoords='data',
            arrowprops=dict(arrowstyle="simple, head_width=0.5, head_length=1, tail_width=0.0001",
                            connectionstyle="arc3", ec ="g", fc = "g", ),
            )

axes.grid(True)
figure.savefig("d_g_versus_d_a.pdf")


