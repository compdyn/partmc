#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import scipy.io, numpy


ncf = scipy.io.netcdf_file("out/urban_plume_process.nc")
time = ncf.variables["time"].data / 3600

num_frac_wet = ncf.variables["num_frac_wet"].data.copy()

print(num_frac_wet)

num_frac_wet[0] = 1.

(figure, axes) = mpl_helper.make_fig(right_margin=0.5)
axes.plot(time, num_frac_wet*100., "g-", label="base")
axes.set_xlabel(r"time / h")
axes.set_ylabel(r"fraction of wet particles / \%")
axes.grid(True)
#axes.legend(loc=(1.05,0))
figure.savefig("figs/num_frac_wet.pdf")

