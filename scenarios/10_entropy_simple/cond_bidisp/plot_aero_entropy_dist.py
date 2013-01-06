#!/usr/bin/env python

import sys, os
sys.path.append("../../../tool")
import mpl_helper
import partmc
import scipy.io
import numpy as np

dist_array = np.zeros([25,200])
i_counter = 0
for (filename, index) in partmc.get_filename_list('out/', r'urban_plume_([0-9]+)_process\.nc'):
    print index, i_counter
    (figure, axes) = mpl_helper.make_fig(right_margin=0.8)

    ncf = scipy.io.netcdf_file(filename)
    entropy = ncf.variables["entropy"].data
    entropy_dist = ncf.variables["entropy_dist"].data
    print entropy_dist.shape

    dist_array[i_counter,:] = entropy_dist/entropy_dist.max()
    i_counter += 1

axes.plot(entropy, dist_array[0,:], 'r', label = "t=0 h")
axes.plot(entropy, dist_array[12,:], 'g', label = "t=12 h")
axes.plot(entropy, dist_array[24,:], 'b', label = "t=24 h")

#axes.annotate(r"t = 0", (entropy[7], [7]),
#              verticalalignment="bottom", horizontalalignment="right",
#              bbox = dict(facecolor='white', edgecolor='white'),
#              xytext=(0, 5), textcoords='offset points')
#axes.annotate(r"$\bar{\hat{H}}$", (time[11], entropy_of_avg_part[11]),
#              verticalalignment="bottom", horizontalalignment="right",
#              bbox = dict(facecolor='white', edgecolor='white'),
#              xytext=(0, 5), textcoords='offset points')
#axes.annotate(r"$\gamma$", (time[17], tot_entropy_ratio[17]),
#              verticalalignment="bottom", horizontalalignment="right",
#              bbox = dict(facecolor='white', edgecolor='white'),
#              xytext=(0, 5), textcoords='offset points')

axes.set_xlabel(r"particle mixing entropy $H_i$")
axes.set_ylabel(r"norm. num. conc. $n(H_i)$")
axes.legend()
axes.set_xlim(0,2)
axes.grid(True)

out_filename = "out/urban_plume_entropy_dist.pdf"
figure.savefig(out_filename)
print out_filename
