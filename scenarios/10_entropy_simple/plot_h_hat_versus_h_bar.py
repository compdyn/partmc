#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import scipy.io
import numpy as np
import matplotlib.pyplot
from matplotlib.patches import Polygon
from matplotlib.patches import FancyArrowPatch

(figure, axes) = mpl_helper.make_fig(figure_width=5)

ncf = scipy.io.netcdf_file("emit/out/urban_plume_process.nc")
avg_part_entropy1 = ncf.variables["avg_part_entropy"].data
entropy_of_avg_part1 = ncf.variables["entropy_of_avg_part"].data
ncf.close()

ncf = scipy.io.netcdf_file("emit_mixed/out/urban_plume_process.nc")
avg_part_entropy2 = ncf.variables["avg_part_entropy"].data
entropy_of_avg_part2 = ncf.variables["entropy_of_avg_part"].data
ncf.close()

ncf = scipy.io.netcdf_file("emit_pure/out/urban_plume_process.nc")
avg_part_entropy3 = ncf.variables["avg_part_entropy"].data
entropy_of_avg_part3 = ncf.variables["entropy_of_avg_part"].data
ncf.close()

ncf = scipy.io.netcdf_file("coag/out/urban_plume_process.nc")
avg_part_entropy4 = ncf.variables["avg_part_entropy"].data
entropy_of_avg_part4 = ncf.variables["entropy_of_avg_part"].data
ncf.close()

ncf = scipy.io.netcdf_file("cond_mono/out/urban_plume_process.nc")
avg_part_entropy5 = ncf.variables["avg_part_entropy"].data
entropy_of_avg_part5 = ncf.variables["entropy_of_avg_part"].data
ncf.close()

ncf = scipy.io.netcdf_file("cond_mono_decrease/out/urban_plume_process.nc")
avg_part_entropy6 = ncf.variables["avg_part_entropy"].data
entropy_of_avg_part6 = ncf.variables["entropy_of_avg_part"].data
ncf.close()

ncf = scipy.io.netcdf_file("cond_bidisp/out/urban_plume_process.nc")
avg_part_entropy7 = ncf.variables["avg_part_entropy"].data
entropy_of_avg_part7 = ncf.variables["entropy_of_avg_part"].data
ncf.close()

ncf = scipy.io.netcdf_file("cond_bidisp_decrease/out/urban_plume_process.nc")
avg_part_entropy8 = ncf.variables["avg_part_entropy"].data
entropy_of_avg_part8 = ncf.variables["entropy_of_avg_part"].data
ncf.close()

print avg_part_entropy5[0], avg_part_entropy5[1], entropy_of_avg_part5[0], entropy_of_avg_part5[1]
axes.plot(avg_part_entropy1/1.1, entropy_of_avg_part1/1.1, "b-", linewidth=2)
axes.plot(avg_part_entropy2/1.4, entropy_of_avg_part2/1.4, "b--", linewidth=2)
axes.plot(avg_part_entropy3/1.1, entropy_of_avg_part3/1.1, "c--", linewidth=2)
axes.plot(avg_part_entropy4/0.7, entropy_of_avg_part4/0.7, "r-", linewidth=2)
axes.plot(avg_part_entropy5/1.4, entropy_of_avg_part5/1.4, "k-", linewidth=2)
axes.plot(avg_part_entropy6/1.4, entropy_of_avg_part6/1.4, "k--", linewidth=2)
axes.plot(avg_part_entropy7/1.4, entropy_of_avg_part7/1.4, "g-", linewidth=2)
axes.plot(avg_part_entropy8/1.4, entropy_of_avg_part8/1.4, "g--", linewidth=2)

axes.set_xlabel(r"$\bar{H}/\ln(A)$")
axes.set_xlim([0,1])
axes.set_xticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
axes.set_ylim([0,1.2])
axes.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1, 1.2])
axes.set_ylabel(r"$\bar{\hat{H}}/\ln{A}$")
pts = np.array([[0,0], [1,1], [0,1]])
p = Polygon(pts,closed=True,alpha=0.05)
axes.add_patch(p)

#axes.annotate("",
#            xy=(avg_part_entropy1[9]/1.1,entropy_of_avg_part1[9]), xycoords='data',
#            xytext=(avg_part_entropy1[5]/1.1,entropy_of_avg_part1[5]), textcoords='data',
#            arrowprops=dict(arrowstyle="simple, head_width=0.5, head_length=1, tail_width=0.001",
#                            connectionstyle="arc3", ec ="b", fc = "b", ),
#            )

#axes.annotate("",
#            xy=(avg_part_entropy2[12]/1.4,entropy_of_avg_part2[12]), xycoords='data',
#            xytext=(avg_part_entropy2[5]/1.4,entropy_of_avg_part2[5]), textcoords='data',
#            arrowprops=dict(arrowstyle="simple, head_width=0.5, head_length=1, tail_width=0.001",
#                            connectionstyle="arc3", ec ="b", fc = "b", ),
#            )

#axes.annotate("",
#            xy=(avg_part_entropy3[15]/1.1,entropy_of_avg_part3[15]), xycoords='data',
#            xytext=(avg_part_entropy3[5]/1.1,entropy_of_avg_part3[5]), textcoords='data',
#            arrowprops=dict(arrowstyle="simple, head_width=0.5, head_length=1, tail_width=0.001",
#                            connectionstyle="arc3", ec ="c", fc = "c", ),
#            )

#axes.annotate("",
#            xy=(avg_part_entropy4[15]/0.7,entropy_of_avg_part4[15]), xycoords='data',
#            xytext=(avg_part_entropy4[5]/0.7,entropy_of_avg_part4[5]), textcoords='data',
#            arrowprops=dict(arrowstyle="simple, head_width=0.5, head_length=1, tail_width=0.001",
#                            connectionstyle="arc3", ec ="r", fc = "r", ),
#            )

#axes.annotate("",
#            xy=(avg_part_entropy5[7]/1.4,entropy_of_avg_part5[7]), xycoords='data',
#            xytext=(avg_part_entropy5[5]/1.4,entropy_of_avg_part5[5]), textcoords='data',
#            arrowprops=dict(arrowstyle="simple, head_width=0.5, head_length=1, tail_width=0.001",
#                            connectionstyle="arc3", ec ="k", fc = "k", ),
#            )

#axes.annotate("",
#            xy=(avg_part_entropy6[17]/1.4,entropy_of_avg_part6[17]), xycoords='data',
#            xytext=(avg_part_entropy6[5]/1.4,entropy_of_avg_part6[5]), textcoords='data',
#            arrowprops=dict(arrowstyle="simple, head_width=0.5, head_length=1, tail_width=0.001",
#                            connectionstyle="arc3", ec ="k", fc = "k", ),
#            )

#axes.annotate("",
#            xy=(avg_part_entropy7[16]/1.4,entropy_of_avg_part7[16]), xycoords='data',
#            xytext=(avg_part_entropy7[6]/1.4,entropy_of_avg_part7[6]), textcoords='data',
#            arrowprops=dict(arrowstyle="simple, head_width=0.5, head_length=1, tail_width=0.001",
#                            connectionstyle="arc3", ec ="g", fc = "g", ),
#            )

#print 'here'
#axes.annotate("",
#            xy=(avg_part_entropy8[17]/1.4,entropy_of_avg_part8[17]), xycoords='data',
#            xytext=(avg_part_entropy8[8]/1.4,entropy_of_avg_part8[8]), textcoords='data',
#            arrowprops=dict(arrowstyle="simple, head_width=0.5, head_length=1, tail_width=0.0001",
#                            connectionstyle="arc3", ec ="g", fc = "g", ),
#            )

axes.grid(True)
figure.savefig("h_hat_versus_h_bar.pdf")
