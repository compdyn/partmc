#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import matplotlib
import scipy.io
import numpy as np
import pandas as pd
 

# Model output
kappa_dist_plot = np.zeros([50,9])

i = 0
for (filename, index) in mpl_helper.get_filename_list('out_pfr_suc2as_2/', r'urban_plume_([0-9]+)_process\.nc'):
    ncf = scipy.io.netcdf_file(filename)
    print(filename)
    num_dist = ncf.variables["num_dist"].data.copy()
    kappa_pmc = ncf.variables["kappa"].data.copy()
    print(kappa_pmc.shape)
    kappa_dist = ncf.variables["kappa_dist"].data.copy() / 1e6
    print(kappa_dist.shape, kappa_pmc.shape)
    kappa_dist_plot[:,i] = kappa_dist
    i = i + 1
    

(figure, axes) = mpl_helper.make_fig(left_margin=1, right_margin=0.5)

axes.plot(kappa_pmc, kappa_dist_plot[:,0]/sum(kappa_dist_plot[:,0]), label="pmc 0 min")
axes.plot(kappa_pmc, kappa_dist_plot[:,1]/sum(kappa_dist_plot[:,1]), label="pmc 10 min")
axes.plot(kappa_pmc, kappa_dist_plot[:,2]/sum(kappa_dist_plot[:,2]), label="pmc 20 min")
#axes.plot(kappa_pmc, kappa_dist_plot[:,3], label="30 min")
axes.plot(kappa_pmc, kappa_dist_plot[:,4]/sum(kappa_dist_plot[:,4]), label="pmc 40 min")
#axes.plot(kappa_pmc, kappa_dist_plot[:,5], label="pmc 50 min")
axes.plot(kappa_pmc, kappa_dist_plot[:,6]/sum(kappa_dist_plot[:,6]), label="pmc 60 min")
#axes.plot(kappa_pmc, kappa_dist_plot[:,7], label="70 min")
#axes.plot(kappa_pmc, kappa_dist_plot[:,8], label="pmc 80 min")

    
axes.set_xscale("linear")
axes.set_xlabel(r"$\kappa$ / 1 ")
axes.set_xlim(0, 0.8)
#axes.set_ylim(5e4, 5e7)

axes.set_yscale("log")
axes.set_ylabel(r"dN/d$\kappa / N_{\rm tot,bin}$ / 1")

axes.grid(True)
axes.legend(loc='upper center')
out_filename = "out_pfr_suc2as_2/urban_plume_kappa_dist.pdf"
figure.savefig(out_filename)
matplotlib.pyplot.close(figure)
