#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import scipy.io, numpy

(figure, axes) = mpl_helper.make_fig(top_margin=0.6,right_margin=0.8)
#axes2 = axes.twinx()

ncf = scipy.io.netcdf_file("out/urban_plume_process.nc")
time = ncf.variables["time"].data.copy() / 3600
num_conc = ncf.variables["tot_num_conc"].data.copy() / 1e6
num_conc_err = ncf.variables["tot_num_conc_ci_offset"].data.copy() / 1e6
mass_conc = ncf.variables["tot_mass_conc"].data.copy() * 1e9
mass_conc_err = ncf.variables["tot_mass_conc_ci_offset"].data.copy() * 1e9

axes.plot(time, num_conc, "b-", label='number conc.')
#axes.errorbar(time, num_conc, num_conc_err, fmt="b-", label='number conc.')
#axes2.errorbar(time, mass_conc, mass_conc_err, fmt="r-", label='mass conc.')
axes.set_xlabel(r"time / h")
axes.set_ylabel(r"num. conc. / $\rm cm^{-3}$")
#axes2.set_ylabel(r"mass conc. / $\rm \mu g\ m^{-3}$")
axes.legend(loc='lower center', bbox_to_anchor=(0.2, 1.0))
#axes2.legend(loc='lower center', bbox_to_anchor=(0.6, 1.0))
axes.grid(True)

out_filename = "out/urban_plume_total.pdf"
print("Writing %s" % out_filename)
figure.savefig(out_filename)
