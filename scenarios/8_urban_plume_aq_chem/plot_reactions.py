#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import scipy.io
import numpy as np

ncf = scipy.io.netcdf_file("out/urban_plume_aq_chem_b_process.nc")
time = ncf.variables["time"].data / 60

forw = []
back = []
conv = []
for i in range(1,12):
    ncf = scipy.io.netcdf_file(f"out/urban_plume_aq_chem_b_0001_{i:08d}.nc")
    forw.append(ncf.variables["aq_chem_rates_forward"].data)
    back.append(ncf.variables["aq_chem_rates_backward"].data)
    conv.append(ncf.variables["aq_chem_rates_conv_factors"].data)
    ncf.close()


forw = np.array(forw) # time x reaction x particle
back = np.array(back)
conv = np.array(conv)

print(forw.shape)

for i_part in range(forw.shape[2]):
    (figure, axes) = mpl_helper.make_fig(right_margin=1.2)
    for i_spec in [0, 1, 2, 4, 9, 12]:
        axes.semilogy(time, forw[:,i_spec,i_part], label=i_spec)

    axes.set_xlabel(r"time / min")
    axes.set_ylabel(r"reaction rate")
    axes.grid(True)
    axes.legend(loc=(1.05,0))
    figure.savefig(f"figs/urban_plume_forw_{i_part}.pdf")

for i_part in range(back.shape[2]):
    (figure, axes) = mpl_helper.make_fig(right_margin=1.2)
    for i_spec in [0, 1, 2, 4, 9, 12]:
        axes.semilogy(time, back[:,i_spec,i_part], label=i_spec)

    axes.set_xlabel(r"time / min")
    axes.set_ylabel(r"reaction rate")
    axes.grid(True)
    axes.legend(loc=(1.05,0))
    figure.savefig(f"figs/urban_plume_back_{i_part}.pdf")

    
