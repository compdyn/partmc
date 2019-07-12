#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import scipy.io
import numpy as np

ncf = scipy.io.netcdf_file("out/urban_plume_aq_chem_b_process.nc")
time = ncf.variables["time"].data / 60

ncf = scipy.io.netcdf_file(f"out/urban_plume_aq_chem_b_0001_00000001.nc")
aero_species = str(ncf.variables["aero_species"].names).split(",")
print(", ".join([f'{i_spec}:{aero_species[i_spec]}' for i_spec in range(len(aero_species))]))
ncf.close()

masses = []
for i in range(1,12):
    ncf = scipy.io.netcdf_file(f"out/urban_plume_aq_chem_b_0001_{i:08d}.nc")
    masses.append(ncf.variables["aero_particle_mass"].data)
    ncf.close()

masses = np.array(masses) # time x species x particle

for i_part in range(masses.shape[2]):
    (figure, axes) = mpl_helper.make_fig(right_margin=1.2)
    for i_spec in [0, 1, 2, 4, 9, 12]:
        axes.semilogy(time, masses[:,i_spec,i_part], label=aero_species[i_spec])

    axes.set_xlabel(r"time / min")
    axes.set_ylabel(r"kg")
    axes.grid(True)
    axes.legend(loc=(1.05,0))
    figure.savefig(f"figs/urban_plume_aq_chem_particle_{i_part}.pdf")
