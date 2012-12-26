#!/usr/bin/env python

import sys, os
sys.path.append("../../../tool")
import mpl_helper
import scipy.io, numpy

(figure, axes) = mpl_helper.make_fig(right_margin=0.8)
axes2 = axes.twinx()

ncf = scipy.io.netcdf_file("out/urban_plume_process.nc")
time = ncf.variables["time"].data / 3600
num_conc = ncf.variables["tot_num_conc"].data / 1e6
num_conc_err = ncf.variables["tot_num_conc_ci_offset"].data / 1e6
tot_entropy_conc = ncf.variables["tot_entropy_conc"].data
tot_entropy_conc_err = ncf.variables["tot_entropy_conc_ci_offset"].data
tot_entropy_of_avg_conc = ncf.variables["tot_entropy_of_avg_conc"].data
tot_entropy_of_avg_conc_err = ncf.variables["tot_entropy_of_avg_conc_ci_offset"].data
tot_entropy_ratio = ncf.variables["tot_entropy_ratio"].data

#axes.errorbar(time, num_conc, num_conc_err, fmt="b-")
axes.errorbar(time, tot_entropy_conc, tot_entropy_conc_err, fmt="r-", label = r'$c_{\rm H}$')
axes.errorbar(time, tot_entropy_of_avg_conc, tot_entropy_of_avg_conc_err, fmt="g-", label = r'$c_{\rm \hat{H}}$')
axes2.plot(time, tot_entropy_ratio, "m-", label = r'$R$')
axes.set_xlabel(r"time / h")
axes.set_ylabel(r"entropy conc. / $\rm m^{-3}$")
axes2.set_ylabel(r"entropy ratio")
axes.legend(loc = 'lower right')
axes.grid(True)

figure.savefig("out/urban_plume_entropy_conc.pdf")
