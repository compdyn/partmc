#!/usr/bin/env python

import sys, os
sys.path.append("../../../tool")
import mpl_helper
import partmc
import scipy.io, numpy

for (filename, index) in partmc.get_filename_list('out/', r'urban_plume_([0-9]+)_process\.nc'):
    (figure, axes) = mpl_helper.make_fig(right_margin=0.8)

    ncf = scipy.io.netcdf_file(filename)
    entropy = ncf.variables["entropy"].data
    dist_ratio_to_entropy_of_avg_part = ncf.variables["dist_ratio_to_entropy_of_avg_part"].data
    dist_ratio_to_entropy_of_avg_part_err = ncf.variables["dist_ratio_to_entropy_of_avg_part_ci_offset"].data


    axes.errorbar(entropy, dist_ratio_to_entropy_of_avg_part / dist_ratio_to_entropy_of_avg_part.max(), dist_ratio_to_entropy_of_avg_part_err / dist_ratio_to_entropy_of_avg_part.max(), fmt="b-")
    axes.set_xlabel(r"entropy")
    axes.set_ylabel(r"num. conc.")
    axes.set_ylim(0, 1)
    axes.grid(True)

    out_filename = "out/urban_plume_dist_ratio_hat_%s.pdf" % index
    figure.savefig(out_filename)
    print out_filename
