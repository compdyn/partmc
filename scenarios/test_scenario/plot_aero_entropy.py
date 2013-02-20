#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import scipy.io, numpy

(figure, axes) = mpl_helper.make_fig(right_margin=0.8)
axes2 = axes.twinx()

ncf = scipy.io.netcdf_file("out/urban_plume_process.nc")
time = ncf.variables["time"].data / 3600
num_conc = ncf.variables["tot_num_conc"].data / 1e6
num_conc_err = ncf.variables["tot_num_conc_ci_offset"].data / 1e6
entropy = ncf.variables["tot_entropy"].data
entropy_err = ncf.variables["tot_entropy_ci_offset"].data
avg_entropy = ncf.variables["tot_entropy_averaged"].data
avg_entropy_err = ncf.variables["tot_entropy_averaged_ci_offset"].data

axes.errorbar(time, num_conc, num_conc_err, fmt="b-")
axes2.errorbar(time, entropy, entropy_err, fmt="r-")
axes2.errorbar(time, avg_entropy, avg_entropy_err, fmt="g-")
axes2.plot(time, entropy / avg_entropy, "m-")
axes.set_xlabel(r"time / h")
axes.set_ylabel(r"num. conc. / $\rm cm^{-3}$")
axes2.set_ylabel(r"entropy")
axes.grid(True)

figure.savefig("out/urban_plume_entropy.pdf")
