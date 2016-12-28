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
d_alpha = ncf.variables["d_alpha"].data
d_alpha_err = ncf.variables["d_alpha_ci_offset"].data
d_gamma = ncf.variables["d_gamma"].data
d_gamma_err = ncf.variables["d_gamma_ci_offset"].data
chi = ncf.variables["chi"].data
chi_err = ncf.variables["chi_ci_offset"].data

axes.errorbar(time, num_conc, num_conc_err, fmt="b-")
axes2.errorbar(time, d_alpha, d_alpha_err, fmt="r-")
axes2.errorbar(time, d_gamma, d_gamma_err, fmt="g-")
axes2.errorbar(time, chi, chi_err, fmt="m-")

axes.set_xlabel(r"time / h")
axes.set_ylabel(r"num. conc. / $\rm cm^{-3}$")
axes2.set_ylabel(r"mixing state metrics")
axes.grid(True)

out_filename = "out/urban_plume_entropy.pdf"
print("Writing %s" % out_filename)
figure.savefig(out_filename)
