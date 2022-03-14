#!/usr/bin/env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import mpl_helper
import scipy.io
import matplotlib

for (filename, index) in mpl_helper.get_filename_list('/Users/nriemer/git/partmc/scenarios/5_tube/out_pfr_suc2as_2/', r'urban_plume_([0-9]+)_process\.nc'):
    (figure, axes) = mpl_helper.make_fig(left_margin=0.7, right_margin=1)

    ncf = scipy.io.netcdf_file(filename, mmap=False)
    sc_2d_dist   = ncf.variables["diam_sc_dist"].data * 1e-6
    wet_num_dist = ncf.variables["num_dist"].data * 1e-6
    sc_grid      = ncf.variables["sc"].data * 100 # Supersaturation level in %
    diam         = ncf.variables["diam"].data * 1e6 # diameter in um
    print(sc_grid)

    
    print(sc_grid[25])
    print(sc_grid[37])
    print(sc_grid[48])
    
    # Calculate the cumulative sum of each ss level
    ccn = np.nancumsum(sc_2d_dist, axis=0)*np.log(sc_grid[1]/sc_grid[0]) 

    for i in (25,37,48):
        ss = sc_grid[i]
        if (i == 25):
            axes.plot(diam, ccn[i]/wet_num_dist, color='g', label=r'$s_{\rm c} = 0.1$\%')
        elif (i == 37):
            axes.plot(diam, ccn[i]/wet_num_dist, color='k', label=r'$s_{\rm c} = 0.3$\%')
        elif (i == 48):
            axes.plot(diam, ccn[i]/wet_num_dist, color='b', label=r'$s_{\rm c} = 0.86$\%')

        axes.legend(loc='lower right')    
        axes.set_xscale('linear')
        axes.set_xlim(1e-2,0.3)
        axes.set_xlabel(r'Diameter, $\rm \mu m$')
        axes.set_ylabel(r'$N_{\rm CCN}/N_{\rm CN}$')
        axes.grid(True)
        out_filename = "/Users/nriemer/git/partmc/scenarios/5_tube/out_pfr_suc2as_2/urban_plume_diam_ccn_dist_lin_%s.pdf" % index
        print(out_filename)
        figure.savefig(out_filename)
        matplotlib.pyplot.close(figure)
        

