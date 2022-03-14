#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import scipy.io, numpy

(figure, axes) = mpl_helper.make_fig(top_margin=0.6,right_margin=0.8, bottom_margin=1.0)
axes2 = axes.twinx()

ncf = scipy.io.netcdf_file("/Users/nriemer/2021/anomitra_project/out_aa_case/urban_plume_process.nc")
out_filename = "/Users/nriemer/2021/anomitra_project/out_aa_case/urban_plume_entropy.pdf"

time = ncf.variables["time"].data.copy() / 3600
d_alpha = ncf.variables["d_alpha"].data.copy()
d_alpha_err = ncf.variables["d_alpha_ci_offset"].data.copy()
d_gamma = ncf.variables["d_gamma"].data.copy()
d_gamma_err = ncf.variables["d_gamma_ci_offset"].data.copy()
chi = ncf.variables["chi"].data.copy()
chi_err = ncf.variables["chi_ci_offset"].data.copy()

print(chi)
print(d_alpha)
print(d_gamma)

axes.plot(time, chi*100, "b-", label=r'$\chi$')
axes2.plot(time, d_alpha, "r-", label=r'$D_{\alpha}$')
axes2.plot(time, d_gamma, "g-", label=r'$D_{\gamma}$')

axes.set_xlabel(r"time / h")
axes.set_ylabel(r"mixing state index in \%")
axes2.set_ylabel(r"mixing state metrics / eff. numb. species")
axes.legend(loc='lower center', bbox_to_anchor=(0.2, -0.5))
axes2.legend(loc='lower center', ncol=2, bbox_to_anchor=(0.6, -0.5))
axes.grid(True)
axes.set_title("Adipic acid to AS")


print("Writing %s" % out_filename)
figure.savefig(out_filename)
