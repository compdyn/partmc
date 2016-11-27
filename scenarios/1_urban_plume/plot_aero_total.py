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
mass_conc = ncf.variables["tot_mass_conc"].data * 1e9
mass_conc_err = ncf.variables["tot_mass_conc_ci_offset"].data * 1e9

axes.errorbar(time, num_conc, num_conc_err, fmt="b-")
axes2.errorbar(time, mass_conc, mass_conc_err, fmt="r-")
axes.set_xlabel(r"time / h")
axes.set_ylabel(r"num. conc. / $\rm cm^{-3}$")
axes2.set_ylabel(r"mass conc. / $\rm \mu g\ m^{-3}$")
axes.grid(True)

out_filename = "out/urban_plume_total.pdf"
print("Writing %s" % out_filename)
figure.savefig(out_filename)
