#!/usr/bin/env python

import sys, os
sys.path.append("../../tool")
import mpl_helper
import scipy.io
import numpy as np

ncf = scipy.io.netcdf_file("out/out_mono_nh32_so22/urban_plume_aq_chem_mono_process.nc")
time = ncf.variables["time"].data / 60

ncf = scipy.io.netcdf_file(f"out/out_mono_nh32_so22/urban_plume_aq_chem_mono_0001_00000001.nc")
aero_species = str(ncf.variables["aero_species"].names).split(",")
print(", ".join([f'{i_spec}:{aero_species[i_spec]}' for i_spec in range(len(aero_species))]))
ncf.close()

masses = []
for i in range(1,602):
    ncf = scipy.io.netcdf_file(f"out/out_mono_nh32_so22/urban_plume_aq_chem_mono_0001_{i:08d}.nc")
    masses.append(ncf.variables["aero_particle_mass"].data.copy())
    ncf.close()

masses1 = np.array(masses) # time x species x particle

Hp_masses1 = masses1[:,60,:]
h2o_masses1 =  masses1[:,19,:]
ph_particle1 = -np.log10(Hp_masses1 * 1000 / h2o_masses1)

ncf = scipy.io.netcdf_file("out/out_mono_oh/urban_plume_aq_chem_mono_process.nc")
time = ncf.variables["time"].data / 60

ncf = scipy.io.netcdf_file(f"out/out_mono_oh/urban_plume_aq_chem_mono_0001_00000001.nc")
aero_species = str(ncf.variables["aero_species"].names).split(",")
ncf.close()

masses = []
for i in range(1,602):
    ncf = scipy.io.netcdf_file(f"out/out_mono_oh/urban_plume_aq_chem_mono_0001_{i:08d}.nc")
    masses.append(ncf.variables["aero_particle_mass"].data.copy())
    ncf.close()

masses2 = np.array(masses) # time x species x particle

Hp_masses2 = masses2[:,60,:]
h2o_masses2 =  masses2[:,19,:]
ph_particle2 = -np.log10(Hp_masses2 * 1000 / h2o_masses2)

for i_part in range(1):#(masses.shape[2]):
    (figure, axes) = mpl_helper.make_fig(right_margin=1.5)
    for i_spec in [23]:
        axes.plot(time, masses1[:,i_spec,i_part],"-", label=aero_species[i_spec]+" nh32so22")
        axes.plot(time, masses2[:,i_spec,i_part],"--", label=aero_species[i_spec]+" oh")
    axes.set_xlabel(r"time / min")
    axes.set_ylabel(r"kg")
    axes.set_xlim([0,10])
    axes.grid(True)
    axes.legend(loc=(1.05,0))
    figure.savefig(f"figs/urban_plume_aq_chem_particle_mono_oh_{i_part}_21_23_29_32.pdf")

for i_part in range(1):#(masses.shape[2]):
    (figure, axes) = mpl_helper.make_fig(right_margin=1.5)
    for i_spec in [63, 78, 0, 3]:
        axes.plot(time, masses1[:,i_spec,i_part], "-", label=aero_species[i_spec]+" nh32so22")
        axes.plot(time, masses2[:,i_spec,i_part], "--", label=aero_species[i_spec]+" oh")
    axes.set_xlabel(r"time / min")
    axes.set_ylabel(r"kg")
    axes.set_xlim([0,10])
    axes.grid(True)
    axes.legend(loc=(1.05,0))
    figure.savefig(f"figs/urban_plume_aq_chem_particle_mono_oh_{i_part}_63_78_0_3.pdf")

for i_part in range(1):#(masses.shape[2]):
    (figure, axes) = mpl_helper.make_fig(right_margin=1.5)
    for i_spec in [19]:
        axes.plot(time, masses1[:,i_spec,i_part], "-", label=aero_species[i_spec]+" nh32so22")
        axes.plot(time, masses2[:,i_spec,i_part], "--", label=aero_species[i_spec]+" oh")
    axes.set_xlabel(r"time / min")
    axes.set_ylabel(r"kg")
    axes.set_xlim([0,10])
    axes.grid(True)
    axes.legend(loc=(1.05,0))
    figure.savefig(f"figs/urban_plume_aq_chem_particle_mono_oh_{i_part}_19.pdf")        

for i_part in range(1):    
    (figure, axes) = mpl_helper.make_fig(right_margin=1.5)
    axes.plot(time, ph_particle1[:,i_part], "-", label='nh32so22')
    axes.plot(time, ph_particle2[:,i_part], "--", label='oh')
    axes.set_xlabel(r"time / min")
    axes.set_ylabel(r"pH")
    axes.set_xlim([0,10])
    axes.grid(True)
    axes.legend(loc=(1.05,0))
    figure.savefig(f"figs/urban_plume_aq_chem_particle_mono_oh_{i_part}_pH.pdf")

    
