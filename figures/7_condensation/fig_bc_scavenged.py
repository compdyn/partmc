#!/usr/bin/env python

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc
const = pmc_data_nc.load_constants("../../src/constants.f90")

in_dir = "../../new_cond/out/"
out_filename = "figs/bc_scavenged.pdf" 

time_array = np.linspace(0,48,49)
frac_array = np.zeros((5,len(time_array)))

for counter in range(1,49):
    in_filename = "cond_%02d_ref_0001_00000601.nc" % counter
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_dir+in_filename)
    particles = pmc_data_nc.aero_particle_array_t(ncf)
    ncf.close()

    final_wet_diameter = particles.diameter()
    is_activated = (final_wet_diameter > 3e-6)

    bc = particles.mass(include = ["BC"]) 
    oc = particles.mass(include = ["OC"])
    so4 = particles.mass(include = ["SO4"])
    no3 = particles.mass(include = ["NO3"])
    nh4 = particles.mass(include = ["NH4"])

    total_bc = sum(bc)
    total_oc = sum(oc)
    total_so4 = sum(so4)
    total_no3 = sum(no3)
    total_nh4 = sum(nh4)

    bc_mass_act = sum(bc[is_activated])
    oc_mass_act = sum(oc[is_activated])
    so4_mass_act = sum(so4[is_activated])
    no3_mass_act = sum(no3[is_activated])
    nh4_mass_act = sum(nh4[is_activated])

    frac_array[0,counter] = bc_mass_act / total_bc
    frac_array[1,counter] = oc_mass_act / total_oc
    frac_array[2,counter] = so4_mass_act / total_so4
    frac_array[3,counter] = no3_mass_act / total_no3
    frac_array[4,counter] = nh4_mass_act / total_nh4

    print 'bc ', counter, total_bc, bc_mass_act

plt.clf()
plt.plot(time_array, frac_array[0,:], 'r-', label = 'bc')
plt.plot(time_array, frac_array[1,:], 'b-', label = 'oc')
plt.plot(time_array, frac_array[2,:], 'g-', label = 'so4')
plt.plot(time_array, frac_array[3,:], 'k-', label = 'no3')
plt.plot(time_array, frac_array[4,:], 'm-', label = 'nh4')
plt.legend(loc = 'lower right')
plt.xlabel("time ")
plt.ylabel("scavenged fraction")
fig = plt.gcf()
fig.savefig(out_filename)
