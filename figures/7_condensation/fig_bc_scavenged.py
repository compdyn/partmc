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
frac_array = np.zeros_like(time_array)

for counter in range(1,49):

    in_filename = "cond_%02d_ref_0001_00000601.nc" % counter
    ncf = Scientific.IO.NetCDF.NetCDFFile(in_dir+in_filename)
    particles = pmc_data_nc.aero_particle_array_t(ncf)
    ncf.close()

    final_wet_diameter = particles.diameter()
    is_activated = (final_wet_diameter > 3e-6)

    bc = particles.mass(include = ["BC"]) 
    total_bc = sum(bc)
    bc_mass_act = sum(bc[is_activated])
    frac_array[counter] = bc_mass_act / total_bc
    print 'bc ', counter, total_bc, bc_mass_act

plt.clf()
plt.plot(time_array, frac_array, 'r-')
plt.xlabel("time ")
plt.ylabel("scavenged BC fraction")
fig = plt.gcf()
fig.savefig(out_filename)
