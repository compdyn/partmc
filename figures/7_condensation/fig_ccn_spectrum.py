#!/usr/bin/env python2.5

import Scientific.IO.NetCDF
import sys
import numpy as np
import matplotlib
matplotlib.use("PDF")
import matplotlib.pyplot as plt
sys.path.append("../../tool")
import pmc_data_nc
const = pmc_data_nc.load_constants("../../src/constants.f90")

in_filename = "../../new_cond/start/urban_plume_wc_0001_00000002.nc"
out_filename = "figs/ccn_spectrum.pdf"
ncf = Scientific.IO.NetCDF.NetCDFFile(in_filename)
particles = pmc_data_nc.aero_particle_array_t(ncf)
env_state = pmc_data_nc.env_state_t(ncf)
ncf.close()

s_crit = (particles.kappa_rh(env_state, const) - 1)*100
total_number = len(particles.id)
ss_array = np.logspace(-2,1,100)
frac_array = np.zeros_like(ss_array)

for i in range(len(ss_array)):
    activated = (s_crit < ss_array[i])
    number_act = sum(activated)
    frac_array[i] = float(number_act) / total_number

plt.grid(True)
plt.semilogx(ss_array, frac_array)
plt.ylabel("CCN/CN")
plt.xlabel("critical supersaturation (%)")
fig = plt.gcf()
fig.savefig(out_filename)

