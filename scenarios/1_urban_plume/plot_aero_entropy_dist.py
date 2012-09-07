#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import partmc
import scipy.io, numpy

for (filename, index) in partmc.get_filename_list('out/', r'urban_plume_([0-9]+)_process\.nc'):
    (figure, axes) = mpl_helper.make_fig(right_margin=0.8)

    ncf = scipy.io.netcdf_file(filename)
    entropy = ncf.variables["entropy"].data
    entropy_dist = ncf.variables["entropy_dist"].data
    entropy_err = ncf.variables["entropy_dist_ci_offset"].data

    entropy_dist_ratio = ncf.variables["entropy_dist_ratio"].data
    entropy_ratio_err = ncf.variables["entropy_dist_ratio_ci_offset"].data

#    axes.errorbar(entropy, entropy_dist / entropy_dist.max(), entropy_err / entropy_dist.max(), fmt="b-")

    axes.errorbar(entropy, entropy_dist_ratio / entropy_dist_ratio.max(), entropy_ratio_err / entropy_dist_ratio.max(), fmt="b-")
    axes.set_xlabel(r"entropy ratio")
    axes.set_ylabel(r"norm. num. conc.")
    axes.set_ylim(0, 1)
    axes.grid(True)

    out_filename = "out/urban_plume_entropy_ratio_%s.pdf" % index
    figure.savefig(out_filename)
    print out_filename
