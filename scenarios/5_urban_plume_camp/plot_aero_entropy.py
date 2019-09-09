#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import scipy.io, numpy

(figure, axes) = mpl_helper.make_fig(top_margin=0.6,right_margin=0.8)
axes2 = axes.twinx()

ncf = scipy.io.netcdf_file("out/urban_plume_process.nc")
time = ncf.variables["time"].data.copy() / 3600
d_alpha = ncf.variables["d_alpha"].data.copy()
d_alpha_err = ncf.variables["d_alpha_ci_offset"].data.copy()
d_gamma = ncf.variables["d_gamma"].data.copy()
d_gamma_err = ncf.variables["d_gamma_ci_offset"].data.copy()
chi = ncf.variables["chi"].data.copy()
chi_err = ncf.variables["chi_ci_offset"].data.copy()

axes.errorbar(time, chi, chi_err, fmt="b-", label=r'$\chi$')
axes2.errorbar(time, d_alpha, d_alpha_err, fmt="r-", label=r'$D_{\alpha}$')
axes2.errorbar(time, d_gamma, d_gamma_err, fmt="g-", label=r'$D_{\gamma}$')

axes.set_xlabel(r"time / h")
axes.set_ylabel(r"mixing state index")
axes2.set_ylabel(r"mixing state metrics / eff. numb. species")
axes.legend(loc='lower center', bbox_to_anchor=(0.2, 1.0))
axes2.legend(loc='lower center', ncol=2, bbox_to_anchor=(0.6, 1.0))
axes.grid(True)

out_filename = "out/urban_plume_entropy.pdf"
print("Writing %s" % out_filename)
figure.savefig(out_filename)
