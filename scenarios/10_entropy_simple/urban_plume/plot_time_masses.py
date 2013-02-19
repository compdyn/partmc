#!/usr/bin/env python

import sys, os
sys.path.append("../../../tool")
import mpl_helper
import scipy.io, numpy

(figure, axes_array) = mpl_helper.make_fig_array(2,1,figure_width=4,vert_sep=0.3,)

ncf = scipy.io.netcdf_file("out/urban_plume2_process.nc")
time = ncf.variables["time"].data / 3600
mass_conc = ncf.variables["tot_mass_conc"].data *1e9
so4_conc = ncf.variables["tot_so4_conc"].data *1e9
no3_conc = ncf.variables["tot_no3_conc"].data *1e9
nh4_conc = ncf.variables["tot_nh4_conc"].data *1e9
bc_conc = ncf.variables["tot_bc_conc"].data *1e9
oc_conc= ncf.variables["tot_oc_conc"].data *1e9
soa_conc = ncf.variables["tot_soa_conc"].data *1e9

axes = axes_array[0][0]
axes.plot(time, so4_conc, "r-")
axes.plot(time, no3_conc, "g-")
axes.plot(time, nh4_conc, "k-")


axes.annotate(r"$\rm SO_4$", (time[202], so4_conc[202]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$\rm NO_3$", (time[50], no3_conc[50]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"$\rm NH_4$", (time[150], nh4_conc[150]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')

axes.set_xlabel(r"time / h")
axes.set_xlim([0,48])
axes.set_xticks([0, 12, 24, 36, 48])
axes.set_ylabel(r"mass conc. / $\rm \mu g \, m^{-3}$")
axes.grid(True)

axes = axes_array[1][0]
axes.plot(time, bc_conc, "r-")
axes.plot(time, oc_conc, "g-")
axes.plot(time, soa_conc,"k-")

axes.annotate(r"BC", (time[74], bc_conc[74]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"OC", (time[180], oc_conc[180]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')
axes.annotate(r"SOA", (time[162], soa_conc[162]),
              verticalalignment="bottom", horizontalalignment="right",
              bbox = dict(facecolor='white', edgecolor='white'),
              xytext=(0, 5), textcoords='offset points')

axes.set_ylabel(r"mass conc. / $\rm \mu g \, m^{-3}$")
axes.grid(True)

mpl_helper.remove_fig_array_axes(axes_array)
figure.savefig("out/urban_plume_time_masses.pdf")
