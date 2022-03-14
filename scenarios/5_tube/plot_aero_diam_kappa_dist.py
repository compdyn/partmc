#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import matplotlib
import scipy.io
import numpy as np

for (filename, index) in mpl_helper.get_filename_list('/Users/nriemer/git/partmc/scenarios/5_tube/out_pfr_suc2as_2/', r'urban_plume_([0-9]+)_process\.nc'):
    (figure, axes, cbar_axes) = mpl_helper.make_fig(left_margin=0.7, right_margin=1, colorbar=True)

    ncf = scipy.io.netcdf_file(filename)
    diam_edges = ncf.variables["diam_edges"].data.copy() * 1e6
    kappa_edges = ncf.variables["kappa_edges"].data.copy()
    diam_kappa_dist = ncf.variables["diam_kappa_dist"].data.copy() * 1e-6
    sc_grid      = ncf.variables["sc"].data * 100 # Supersaturation level in %

    print("sc_grid")
    print(sc_grid)

    sc1 = 1+0.3e-2
    sc2 = 1+0.86e-2
    sc3 = 1+0.1e-2
    sigma = 0.07275
    M_w   = 18e-3
    rho_w = 1000
    T     = 293
    R     = 8.3145
    A     = 4e0 * sigma * M_w / R / T / rho_w
    sc_line1 = 4e0 * A**3e0 / (27e0 * (diam_edges / 1e6)**3e0 * (np.log(sc1))**2e0)
    sc_line2 = 4e0 * A**3e0 / (27e0 * (diam_edges / 1e6)**3e0 * (np.log(sc2))**2e0)
    sc_line3 = 4e0 * A**3e0 / (27e0 * (diam_edges / 1e6)**3e0 * (np.log(sc3))**2e0)

    p = axes.pcolor(diam_edges, kappa_edges, diam_kappa_dist,
                    norm = matplotlib.colors.LogNorm(vmin=1e3, vmax=1e6), linewidths = 0.1)
    axes.plot(diam_edges, sc_line1, 'k')
    axes.plot(diam_edges, sc_line2, 'b')
    axes.plot(diam_edges, sc_line3, 'g')
    axes.set_xscale("linear")
    axes.set_xlabel(r"dry diameter $D_{\rm dry}$ / $\rm \mu m$")
    axes.set_xlim(1e-2, 0.3)

    axes.set_yscale("linear")
    axes.set_ylabel(r"hygroscopicity parameter $\kappa$ / 1")
    axes.set_ylim(0, 1)

    axes.annotate(r'0.3\%',  xy=(0.07, 1), xycoords='data', 
            xytext=(0.06, 1), textcoords='data',
            horizontalalignment='center', verticalalignment='bottom',
            )

    axes.annotate(r'0.86\%',  xy=(0.03, 1), xycoords='data', 
            xytext=(0.03, 1), textcoords='data',
            horizontalalignment='center', verticalalignment='bottom',
            )

    axes.annotate(r'0.1\%',  xy=(0.12, 1), xycoords='data', 
            xytext=(0.12, 1), textcoords='data',
            horizontalalignment='center', verticalalignment='bottom',
            )

    axes.grid(True)
    cbar = figure.colorbar(p, cax=cbar_axes, format=matplotlib.ticker.LogFormatterMathtext(),
                           orientation='vertical')
    cbar_axes.xaxis.set_label_position('top')
    cbar.set_label(r"number conc. $n(D_{\rm dry},\kappa)$ / $\rm cm^{-3}$")

    out_filename = "/Users/nriemer/git/partmc/scenarios/5_tube/out_pfr_suc2as_2/urban_plume_diam_kappa_dist_lin_%s.pdf" % index
    print("Writing %s" % out_filename)
    figure.savefig(out_filename)
    matplotlib.pyplot.close(figure)
