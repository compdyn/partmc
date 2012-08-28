#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import scipy.io, numpy

(figure, axes) = mpl_helper.make_fig(right_margin=0.8)

ncf = scipy.io.netcdf_file("out/urban_plume_comp_process.nc")
time = ncf.variables["time"].data / 3600
ccn01 = ncf.variables["ccn_01"].data / 1e6
ccn01_err = ncf.variables["ccn_01_ci_offset"].data / 1e6
ccn03 = ncf.variables["ccn_03"].data / 1e6
ccn03_err = ncf.variables["ccn_03_ci_offset"].data / 1e6
ccn06 = ncf.variables["ccn_06"].data / 1e6
ccn06_err = ncf.variables["ccn_06_ci_offset"].data / 1e6

axes.errorbar(time, ccn01, ccn01_err, fmt="b-", label = '0.1 \%')
axes.errorbar(time, ccn03, ccn03_err, fmt="g-", label = '0.3 \%')
axes.errorbar(time, ccn06, ccn06_err, fmt="r-", label = '0.6 \%')
axes.legend()
axes.set_xlabel(r"time / h")
axes.set_ylabel(r"CCN conc. / $\rm cm^{-3}$")
axes.grid(True)

figure.savefig("out/urban_plume_ccn_comp.pdf")
