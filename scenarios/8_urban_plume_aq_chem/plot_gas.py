#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import scipy.io
import numpy as np

ncf = scipy.io.netcdf_file("out/out_mono_nh32/urban_plume_aq_chem_mono_process.nc")
time = ncf.variables["time"].data / 60

ncf = scipy.io.netcdf_file(f"out/out_mono_nh32/urban_plume_aq_chem_mono_0001_00000001.nc")
gas_species = str(ncf.variables["gas_species"].names).split(",")
print(", ".join([f'{i_spec}:{gas_species[i_spec]}' for i_spec in range(len(gas_species))]))
ncf.close()

gas = []
for i in range(1,602):
    ncf = scipy.io.netcdf_file(f"out/out_mono_nh32/urban_plume_aq_chem_mono_0001_{i:08d}.nc")
    gas.append(ncf.variables["gas_mixing_ratio"].data.copy())
    ncf.close()

gas1 = np.array(gas) # time x species

ncf = scipy.io.netcdf_file("out/out_mono_nh32_so22/urban_plume_aq_chem_mono_process.nc")
time = ncf.variables["time"].data / 60

ncf = scipy.io.netcdf_file(f"out/out_mono_nh32_so22/urban_plume_aq_chem_mono_0001_00000001.nc")
gas_species = str(ncf.variables["gas_species"].names).split(",")
ncf.close()

gas = []
for i in range(1,602):
    ncf = scipy.io.netcdf_file(f"out/out_mono_nh32_so22/urban_plume_aq_chem_mono_0001_{i:08d}.nc")
    gas.append(ncf.variables["gas_mixing_ratio"].data.copy())
    ncf.close()

gas2 = np.array(gas) # time x species


(figure, axes) = mpl_helper.make_fig(right_margin=2)
for i_spec in [3, 17]:
    axes.plot(time, gas1[:,i_spec], label=gas_species[i_spec]+" nh32")
    axes.plot(time, gas2[:,i_spec], label=gas_species[i_spec]+" nhnh32so22")
    
    axes.set_xlabel(r"time / min")
    axes.set_ylabel(r"mixing ratio / ppb")
    axes.grid(True)
    axes.legend(loc=(1.05,0))
    figure.savefig(f"figs/urban_plume_mono_nh32_so22_gas.pdf")
