#!/usr/bin/env python

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import partmc
const = partmc.constants_t("../../src/constants.f90")

in_filename = "../../scenarios/1_urban_plume/out/urban_plume_wc_0001_00000025.nc"
out_filename = "figs/ccn_spectrum_bc_wc_24.pdf"
ncf = Scientific.IO.NetCDF.NetCDFFile(in_filename)
particles = partmc.aero_particle_array_t(ncf)
env_state = partmc.env_state_t(ncf)
ncf.close()

s_crit = (particles.critical_rel_humids(env_state, const) - 1)*100
bc = particles.masses(include = ["BC"])
total_bc = sum(bc)
ss_array = np.logspace(-2,1,100)
frac_array = np.zeros_like(ss_array)

for i in range(len(ss_array)):
    activated = (s_crit < ss_array[i])
    bc_mass_act = sum(bc[activated])
    frac_array[i] = bc_mass_act / total_bc

plt.grid(True)
plt.semilogx(ss_array, frac_array)
plt.ylabel("scavenged BC fraction")
plt.xlabel("critical supersaturation (%)")
fig = plt.gcf()
fig.savefig(out_filename)

