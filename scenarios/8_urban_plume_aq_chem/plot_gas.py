#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import scipy.io
import numpy as np

ncf = scipy.io.netcdf_file("out/urban_plume_aq_chem_b_process.nc")
time = ncf.variables["time"].data / 60

ncf = scipy.io.netcdf_file(f"out/urban_plume_aq_chem_b_0001_00000001.nc")
gas_species = str(ncf.variables["gas_species"].names).split(",")
print(", ".join([f'{i_spec}:{gas_species[i_spec]}' for i_spec in range(len(gas_species))]))
ncf.close()

gas = []
for i in range(1,12):
    ncf = scipy.io.netcdf_file(f"out/urban_plume_aq_chem_b_0001_{i:08d}.nc")
    gas.append(ncf.variables["gas_mixing_ratio"].data)
    ncf.close()


gas = np.array(gas) # time x species
print(gas.shape)    

for i_part in range(gas.shape[1]):
    (figure, axes) = mpl_helper.make_fig(right_margin=1.2)
    for i_spec in [0, 1, 2, 4, 9, 12]:
        axes.semilogy(time, gas[:,i_spec], label=gas_species[i_spec])

    axes.set_xlabel(r"time / min")
    axes.set_ylabel(r"ppb")
    axes.grid(True)
    axes.legend(loc=(1.05,0))
    figure.savefig(f"figs/urban_plume_gas.pdf")
