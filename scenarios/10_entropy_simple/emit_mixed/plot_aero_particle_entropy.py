#!/usr/bin/env python

import sys, os
sys.path.append("../../../tool")
import mpl_helper
import scipy.io, numpy

(figure, axes) = mpl_helper.make_fig(right_margin=0.8)
axes2 = axes.twinx()

ncf = scipy.io.netcdf_file("out/urban_plume_process.nc")
time = ncf.variables["time"].data / 3600

avg_part_entropy = ncf.variables["avg_part_entropy"].data
avg_part_entropy_err = ncf.variables["avg_part_entropy_ci_offset"].data
entropy_of_avg_part = ncf.variables["entropy_of_avg_part"].data
entropy_of_avg_part_err = ncf.variables["entropy_of_avg_part_ci_offset"].data
tot_entropy_ratio = ncf.variables["tot_entropy_ratio"].data

axes.errorbar(time, avg_part_entropy, avg_part_entropy_err, fmt="r-", label = r'$\bar{H}$')
axes.errorbar(time, entropy_of_avg_part, entropy_of_avg_part_err, fmt="g-", label = r'$\bar{\hat{H}}$')
axes2.plot(time, tot_entropy_ratio, "m-", label = r'$R$')
axes.set_xlabel(r"time / h")
axes.set_ylabel(r"entropy")
axes2.set_ylabel(r"entropy ratio")
axes.legend(loc = 'lower right')
axes.grid(True)

figure.savefig("out/urban_plume_particle_entropy.pdf")
