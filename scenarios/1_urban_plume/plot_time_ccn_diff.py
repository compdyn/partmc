#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import scipy.io, numpy

(figure, axes) = mpl_helper.make_fig(top_margin=0.4,right_margin=0.8)

ncf1 = scipy.io.netcdf_file("out/urban_plume_process.nc")
time1 = ncf1.variables["time"].data / 3600
entropy1 = ncf1.variables["tot_entropy"].data
entropy_err1 = ncf1.variables["tot_entropy_ci_offset"].data
avg_entropy1 = ncf1.variables["tot_entropy_averaged"].data
avg_entropy_err1 = ncf1.variables["tot_entropy_averaged_ci_offset"].data

ncf2 = scipy.io.netcdf_file("out/urban_plume_comp_process.nc")
time2 = ncf2.variables["time"].data / 3600
entropy2 = ncf2.variables["tot_entropy"].data
entropy_err2 = ncf2.variables["tot_entropy_ci_offset"].data
avg_entropy2 = ncf2.variables["tot_entropy_averaged"].data
avg_entropy_err2 = ncf2.variables["tot_entropy_averaged_ci_offset"].data


ccn01_1 = ncf1.variables["ccn_01"].data / 1e6
ccn03_1 = ncf1.variables["ccn_03"].data / 1e6
ccn06_1 = ncf1.variables["ccn_06"].data / 1e6

ccn01_2 = ncf2.variables["ccn_01"].data / 1e6
ccn03_2 = ncf2.variables["ccn_03"].data / 1e6
ccn06_2 = ncf2.variables["ccn_06"].data / 1e6

label = "Difference (comp(25 bin) - ref)/ ref"
axes.plot(time1, (ccn01_2-ccn01_1)/ccn01_1, "b-")
axes.plot(time1, (ccn03_2-ccn03_1)/ccn03_1, "g-")
axes.plot(time1, (ccn06_2-ccn06_1)/ccn06_1, "r-")
axes.set_xlabel(r"time / h")
axes.set_ylabel(r"rel. CCN conc. diff. / 1 ")
axes.set_title(label)
axes.grid(True)

figure.savefig("out/urban_plume_time_ccn_diff.pdf")


